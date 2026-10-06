// SABDB vs HBSAE vs direct estimation: in-browser simulation.
// Gibbs samplers mirror the priors in github.com/Sinan-Yavuz/SABDB:
//   SABDB (simulation/SABDB_Sim_unw.stan): beta_s ~ N(mu, diag(tau2)), mu ~ N(0, 3),
//     tau2_k ~ Inv-Gamma(1, 1), sigma_y ~ Half-Cauchy(0, 5); y and covariates standardized.
//   HBSAE (pisa-unit-level/HBSAE_BHF.stan, brms in the simulation): fixed subgroup
//     effects + random state intercept, beta ~ N(0, 10), sigma_u, sigma_y ~ Half-Cauchy(0, 5).
// Half-Cauchy priors use the Makalic & Schmidt (2016) inverse-gamma mixture so every
// update is a conjugate draw.

var SABDBSim = (function () {
  var GROUPS = ["White", "Black", "Hispanic", "Asian"];
  var NAT_MEAN = [292, 260, 268, 310];     // approx. NAEP 2019 grade 8 math means
  var BASE_SHARE = [0.50, 0.16, 0.29, 0.05];
  var WITHIN_SD = 36;
  var G = GROUPS.length, K = G;

  function rng(seed) {
    var t = seed >>> 0;
    var u = function () {
      t += 0x6D2B79F5; var r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
    var spare = null;
    u.normal = function () {
      if (spare !== null) { var s = spare; spare = null; return s; }
      var a, b, q;
      do { a = 2 * u() - 1; b = 2 * u() - 1; q = a * a + b * b; } while (q >= 1 || q === 0);
      var f = Math.sqrt(-2 * Math.log(q) / q); spare = b * f; return a * f;
    };
    u.gamma = function (shape) {
      if (shape < 1) return u.gamma(shape + 1) * Math.pow(u(), 1 / shape);
      var d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
      for (;;) {
        var x, v; do { x = u.normal(); v = 1 + c * x; } while (v <= 0);
        v = v * v * v; var w = u();
        if (w < 1 - 0.0331 * x * x * x * x || Math.log(w) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
      }
    };
    u.invgamma = function (shape, scale) { return scale / u.gamma(shape); };
    return u;
  }

  // Population + sample. Keeps sufficient statistics per state x subgroup cell.
  function simulate(p, r) {
    var S = p.states, cells = [], truth = [];
    for (var s = 0; s < S; s++) {
      var level = r.normal() * p.levelSD, w = [], tot = 0, g;
      for (g = 0; g < G; g++) { w[g] = BASE_SHARE[g] * Math.exp(0.6 * r.normal()); tot += w[g]; }
      var cum = [], c = 0;
      for (g = 0; g < G; g++) { c += w[g] / tot; cum[g] = c; }
      var t = [], cs = [];
      for (g = 0; g < G; g++) {
        t[g] = NAT_MEAN[g] + level + (g === 0 ? 0 : r.normal() * p.gapSD);
        cs[g] = { n: 0, sy: 0, syy: 0 };
      }
      for (var i = 0; i < p.n; i++) {
        var x = r(); g = 0; while (g < G - 1 && x > cum[g]) g++;
        var y = t[g] + WITHIN_SD * r.normal();
        cs[g].n++; cs[g].sy += y; cs[g].syy += y * y;
      }
      truth.push(t); cells.push(cs);
    }
    return { S: S, cells: cells, truth: truth };
  }

  // Standardize y and the dummy columns exactly as the Stan transformed data block does.
  function design(d) {
    var N = 0, sy = 0, syy = 0, ng = new Array(G).fill(0), s, g;
    for (s = 0; s < d.S; s++) for (g = 0; g < G; g++) {
      var c = d.cells[s][g]; N += c.n; sy += c.sy; syy += c.syy; ng[g] += c.n;
    }
    var muY = sy / N, sdY = Math.sqrt((syy - N * muY * muY) / (N - 1));
    var mx = [0], sx = [1];
    for (g = 1; g < G; g++) { var m = ng[g] / N; mx[g] = m; sx[g] = Math.sqrt(m * (1 - m) * N / (N - 1)); }
    var V = [];
    for (g = 0; g < G; g++) {
      var v = [1];
      for (var k = 1; k < K; k++) v[k] = ((g === k ? 1 : 0) - mx[k]) / sx[k];
      V.push(v);
    }
    var st = d.cells.map(function (cs) {
      return cs.map(function (c) {
        return { n: c.n, sy: (c.sy - c.n * muY) / sdY, syy: (c.syy - 2 * muY * c.sy + c.n * muY * muY) / (sdY * sdY) };
      });
    });
    return { N: N, muY: muY, sdY: sdY, sx: sx, V: V, st: st };
  }

  function dot(a, b) { var s = 0; for (var i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

  // Draw from N(A^-1 b, A^-1) for a small symmetric positive-definite A (row-major K x K).
  function mvnCanon(A, b, r) {
    var n = b.length, L = new Float64Array(n * n), i, j, k;
    for (i = 0; i < n; i++) for (j = 0; j <= i; j++) {
      var s = A[i * n + j];
      for (k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = i === j ? Math.sqrt(s) : s / L[j * n + j];
    }
    var z = new Float64Array(n);
    for (i = 0; i < n; i++) { var t = b[i]; for (k = 0; k < i; k++) t -= L[i * n + k] * z[k]; z[i] = t / L[i * n + i]; }
    for (i = 0; i < n; i++) z[i] += r.normal();
    var x = new Array(n);
    for (i = n - 1; i >= 0; i--) { var u = z[i]; for (k = i + 1; k < n; k++) u -= L[k * n + i] * x[k]; x[i] = u / L[i * n + i]; }
    return x;
  }

  function newDraws(S, kept) {
    var a = []; for (var s = 0; s < S; s++) { a.push([]); for (var g = 0; g < G; g++) a[s].push(new Float64Array(kept)); }
    return a;
  }

  function fitSABDB(d, D, r, iters, burn, nu1, nu2) {
    var S = d.S, kept = iters - burn, s, g, k, i, j;
    var beta = [], mu = new Array(K).fill(0), tau2 = new Array(K).fill(1), sig2 = 1, a = 1;
    for (s = 0; s < S; s++) beta.push(new Array(K).fill(0));
    var draws = newDraws(S, kept), tauSum = new Array(K).fill(0);
    for (var it = 0; it < iters; it++) {
      var rss = 0;
      for (s = 0; s < S; s++) {
        var A = new Float64Array(K * K), b = new Array(K);
        for (k = 0; k < K; k++) { A[k * K + k] = 1 / tau2[k]; b[k] = mu[k] / tau2[k]; }
        for (g = 0; g < G; g++) {
          var c = D.st[s][g], v = D.V[g];
          if (!c.n) continue;
          for (i = 0; i < K; i++) { b[i] += v[i] * c.sy / sig2; for (j = 0; j < K; j++) A[i * K + j] += c.n * v[i] * v[j] / sig2; }
        }
        beta[s] = mvnCanon(A, b, r);
        for (g = 0; g < G; g++) {
          var cc = D.st[s][g], m = dot(D.V[g], beta[s]);
          rss += cc.syy - 2 * m * cc.sy + cc.n * m * m;
        }
      }
      sig2 = r.invgamma((D.N + 1) / 2, rss / 2 + 1 / a);
      a = r.invgamma(1, 1 / 25 + 1 / sig2);
      for (k = 0; k < K; k++) {
        var sb = 0; for (s = 0; s < S; s++) sb += beta[s][k];
        var prec = S / tau2[k] + 1 / 9;
        mu[k] = (sb / tau2[k]) / prec + r.normal() / Math.sqrt(prec);
        var ss = 0; for (s = 0; s < S; s++) ss += (beta[s][k] - mu[k]) * (beta[s][k] - mu[k]);
        tau2[k] = r.invgamma(nu1 + S / 2, nu2 + ss / 2);
      }
      if (it >= burn) {
        var t = it - burn;
        for (s = 0; s < S; s++) for (g = 0; g < G; g++) draws[s][g][t] = D.muY + D.sdY * dot(D.V[g], beta[s]);
        for (k = 0; k < K; k++) tauSum[k] += Math.sqrt(tau2[k]) * D.sdY / D.sx[k];
      }
    }
    return { draws: draws, tau: tauSum.map(function (x) { return x / kept; }) };
  }

  function fitHBSAE(d, D, r, iters, burn) {
    var S = d.S, kept = iters - burn, s, g, i, j;
    var beta = new Array(K).fill(0), u = new Array(S).fill(0), sig2 = 1, su2 = 1, a = 1, au = 1;
    var draws = newDraws(S, kept), suSum = 0;
    var XtX = new Float64Array(K * K), ns = new Array(S).fill(0);
    for (s = 0; s < S; s++) for (g = 0; g < G; g++) {
      var c = D.st[s][g], v = D.V[g]; ns[s] += c.n;
      for (i = 0; i < K; i++) for (j = 0; j < K; j++) XtX[i * K + j] += c.n * v[i] * v[j];
    }
    for (var it = 0; it < iters; it++) {
      var A = new Float64Array(K * K), b = new Array(K).fill(0);
      for (i = 0; i < K * K; i++) A[i] = XtX[i] / sig2;
      for (i = 0; i < K; i++) A[i * K + i] += 1 / 100;
      for (s = 0; s < S; s++) for (g = 0; g < G; g++) {
        var cc = D.st[s][g];
        for (i = 0; i < K; i++) b[i] += D.V[g][i] * (cc.sy - cc.n * u[s]) / sig2;
      }
      beta = mvnCanon(A, b, r);
      var fx = D.V.map(function (v) { return dot(v, beta); });
      var su = 0, rss = 0;
      for (s = 0; s < S; s++) {
        var res = 0;
        for (g = 0; g < G; g++) res += D.st[s][g].sy - D.st[s][g].n * fx[g];
        var prec = ns[s] / sig2 + 1 / su2;
        u[s] = (res / sig2) / prec + r.normal() / Math.sqrt(prec);
        su += u[s] * u[s];
        for (g = 0; g < G; g++) {
          var c2 = D.st[s][g], m = fx[g] + u[s];
          rss += c2.syy - 2 * m * c2.sy + c2.n * m * m;
        }
      }
      su2 = r.invgamma((S + 1) / 2, su / 2 + 1 / au);
      au = r.invgamma(1, 1 / 25 + 1 / su2);
      sig2 = r.invgamma((D.N + 1) / 2, rss / 2 + 1 / a);
      a = r.invgamma(1, 1 / 25 + 1 / sig2);
      if (it >= burn) {
        var t = it - burn;
        for (s = 0; s < S; s++) for (g = 0; g < G; g++) draws[s][g][t] = D.muY + D.sdY * (fx[g] + u[s]);
        suSum += Math.sqrt(su2) * D.sdY;
      }
    }
    return { draws: draws, tauLevel: suSum / kept };
  }

  function summarize(arr) {
    var x = Array.prototype.slice.call(arr).sort(function (p, q) { return p - q; });
    var n = x.length, m = 0; for (var i = 0; i < n; i++) m += x[i];
    return { est: m / n, lo: x[Math.floor(0.025 * (n - 1))], hi: x[Math.ceil(0.975 * (n - 1))] };
  }

  // One replication: simulate, fit all three methods, return per-domain results.
  function replicate(p, seed, opt) {
    opt = opt || {};
    var iters = opt.iters || 700, burn = opt.burn || 250;
    var r = rng(seed), d = simulate(p, r), D = design(d);
    var sab = fitSABDB(d, D, r, iters, burn, p.nu1 || 1, p.nu2 || 1), hb = fitHBSAE(d, D, r, iters, burn);
    var rows = [];
    for (var s = 0; s < d.S; s++) for (var g = 0; g < G; g++) {
      var c = d.cells[s][g], dir = null;
      if (c.n >= 2) {
        var m = c.sy / c.n, v = (c.syy - c.n * m * m) / (c.n - 1), se = Math.sqrt(v / c.n);
        dir = { est: m, lo: m - 1.96 * se, hi: m + 1.96 * se };
      }
      rows.push({ s: s, g: g, n: c.n, truth: d.truth[s][g], direct: dir,
        hbsae: summarize(hb.draws[s][g]), sabdb: summarize(sab.draws[s][g]) });
    }
    return { rows: rows, tau: sab.tau, hbTauLevel: hb.tauLevel };
  }

  return { GROUPS: GROUPS, WITHIN_SD: WITHIN_SD, replicate: replicate };
})();

if (typeof module !== "undefined") module.exports = SABDBSim;
