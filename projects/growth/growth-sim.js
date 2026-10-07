// Quadratic growth model: ML vs noninformative Bayes (very wide priors) vs Bayes with priors from an earlier cohort.
// Mirrors the design of Yavuz (2021), Hacettepe University PhD thesis:
//   Level 1: Y_tij = pi0 + pi1 t + pi2 t^2 + e,  t in {0,1,2,3,4,5,7,9,11} (K fall .. grade 5 spring)
//   Level 2: each pi regressed on male and SES; random student intercept and linear slope
//   Level 3: school intercept and slope regressed on public school; random school intercept and slope
// True values come from the thesis's ECLS-K 2011 estimates (tables 5, 14 and 15).
//
// To run in the browser, estimation is two-stage, which is close to the mixed model's fixed
// effects when every student has the same time points (as here):
//   1. per-student OLS growth curve,
//   2. for each growth coefficient, GLS on the student covariates with a school random effect
//      (variance components by method of moments).
// The Bayesian estimates combine that normal likelihood with normal priors (conjugate update).
// Each replication's prior mean is a fresh draw, as if re-estimated from an earlier cohort:
// truth + drift * SD + N(0, SD^2), so a prior with drift 0 is right on average, not exactly right.

var GrowthSim = (function () {
  var TIMES = [0, 1, 2, 3, 4, 5, 7, 9, 11];

  // Fixed effects, rows = [start, linear growth, quadratic], columns = [intercept, male, SES, public].
  // Intercepts are calibrated so the population mean curve matches the thesis's Table 5 means.
  var TRUTH = [
    [38.10, 0.18, 4.45, -2.00],
    [12.43, 0.94, 0.27, -0.03],
    [-0.488, -0.06, -0.01, 0.04]
  ];
  // Prior SDs from the ECLS-K 1998 analysis (thesis Table 11, before inflation).
  var PRIOR_SD = [
    [10, 0.64, 0.50, 1.01],
    [1.83, 0.22, 0.16, 0.27],
    [0.15, 0.02, 0.01, 0.02]
  ];
  // Thesis Table 15 (ML). The school slope variance isn't reported there, so 0.10 is an assumption.
  var VC = { sigma2: 39.23, tauStuInt: 93.01, tauStuSlope: 0.87, corStu: 0.08, tauSchInt: 6.42, tauSchSlope: 0.10 };
  var PUBLIC_SHARE = 153 / 165;   // thesis Table 4

  // Parameters reported in the demo: [equation, column, label]
  var PARAMS = [
    { k: 1, c: 0, key: "t", name: "Growth rate at K entry" },
    { k: 0, c: 2, key: "ses0", name: "SES → starting score" },
    { k: 1, c: 2, key: "ses1", name: "SES → initial growth rate" },
    { k: 0, c: 1, key: "male0", name: "Male → starting score" },
    { k: 1, c: 1, key: "male1", name: "Male → initial growth rate" },
    { k: 0, c: 3, key: "pub0", name: "Public school → starting score" },
    { k: 1, c: 3, key: "pub1", name: "Public school → initial growth rate" }
  ];

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
    return u;
  }

  function chol(A, p) {
    var L = new Float64Array(p * p);
    for (var i = 0; i < p; i++) for (var j = 0; j <= i; j++) {
      var s = A[i * p + j];
      for (var k = 0; k < j; k++) s -= L[i * p + k] * L[j * p + k];
      if (i === j) { if (s <= 1e-12) return null; L[i * p + i] = Math.sqrt(s); }
      else L[i * p + j] = s / L[j * p + j];
    }
    return L;
  }
  function inverse(A, p) {
    var L = chol(A, p); if (!L) return null;
    var inv = new Float64Array(p * p);
    for (var j = 0; j < p; j++) {
      var z = new Float64Array(p), x = new Float64Array(p), i, k;
      for (i = 0; i < p; i++) { var t = i === j ? 1 : 0; for (k = 0; k < i; k++) t -= L[i * p + k] * z[k]; z[i] = t / L[i * p + i]; }
      for (i = p - 1; i >= 0; i--) { var v = z[i]; for (k = i + 1; k < p; k++) v -= L[k * p + i] * x[k]; x[i] = v / L[i * p + i]; }
      for (i = 0; i < p; i++) inv[i * p + j] = x[i];
    }
    return inv;
  }
  function matvec(A, b, p) {
    var out = new Float64Array(p);
    for (var i = 0; i < p; i++) for (var j = 0; j < p; j++) out[i] += A[i * p + j] * b[j];
    return out;
  }

  // (T'T)^-1 T' for the per-student OLS growth curve.
  var PROJ = (function () {
    var T = TIMES.map(function (t) { return [1, t, t * t]; }), TtT = new Float64Array(9);
    T.forEach(function (r) { for (var a = 0; a < 3; a++) for (var b = 0; b < 3; b++) TtT[a * 3 + b] += r[a] * r[b]; });
    var inv = inverse(TtT, 3);
    return T.map(function (r) { return matvec(inv, r, 3); });
  })();

  function simulate(q, r) {
    var J = q.schools, nj = q.perSchool, nPub = Math.min(J - 1, Math.max(1, Math.round(PUBLIC_SHARE * J)));
    var sd = Math.sqrt, students = [];
    var sI = sd(VC.tauStuInt), sS = sd(VC.tauStuSlope), rho = VC.corStu;
    for (var j = 0; j < J; j++) {
      var pub = j < nPub ? 1 : 0, sesSchool = r.normal();
      var u0 = r.normal() * sd(VC.tauSchInt), u1 = r.normal() * sd(VC.tauSchSlope);
      for (var i = 0; i < nj; i++) {
        var male = r() < 0.5 ? 1 : 0, ses = 0.6 * sesSchool + 0.8 * r.normal();
        var x = [1, male, ses, pub], z0 = r.normal(), z1 = rho * z0 + Math.sqrt(1 - rho * rho) * r.normal();
        var pi = TRUTH.map(function (row) { return row[0] * x[0] + row[1] * x[1] + row[2] * x[2] + row[3] * x[3]; });
        pi[0] += u0 + sI * z0; pi[1] += u1 + sS * z1;
        var y = TIMES.map(function (t) { return pi[0] + pi[1] * t + pi[2] * t * t + r.normal() * sd(VC.sigma2); });
        students.push({ school: j, x: x, y: y });
      }
    }
    return { students: students, J: J };
  }

  // GLS with a compound-symmetric school effect, for one growth coefficient.
  function glsEquation(d, k) {
    var p = 4, S = d.students, N = S.length, J = d.J;
    var bk = S.map(function (s) { var v = 0; for (var t = 0; t < TIMES.length; t++) v += PROJ[t][k] * s.y[t]; return v; });
    // OLS to get residuals for the variance components.
    var XtX = new Float64Array(p * p), Xty = new Float64Array(p), a, b, i;
    S.forEach(function (s, i) { for (a = 0; a < p; a++) { Xty[a] += s.x[a] * bk[i]; for (b = 0; b < p; b++) XtX[a * p + b] += s.x[a] * s.x[b]; } });
    var inv0 = inverse(XtX, p); if (!inv0) return null;
    var beta0 = matvec(inv0, Xty, p), res = S.map(function (s, i) { var f = 0; for (a = 0; a < p; a++) f += s.x[a] * beta0[a]; return bk[i] - f; });
    var sum = new Float64Array(J), cnt = new Float64Array(J);
    S.forEach(function (s, i) { sum[s.school] += res[i]; cnt[s.school]++; });
    var within = 0; S.forEach(function (s, i) { var e = res[i] - sum[s.school] / cnt[s.school]; within += e * e; });
    var omega = within / Math.max(1, N - J), grand = 0, between = 0, nbar = N / J;
    for (var j = 0; j < J; j++) grand += sum[j] / cnt[j]; grand /= J;
    for (j = 0; j < J; j++) { var m = sum[j] / cnt[j] - grand; between += m * m; }
    var tau = Math.max(0, between / Math.max(1, J - 1) - omega / nbar);
    // GLS: V_j^-1 = (I - c J) / omega, c = tau / (omega + n_j tau)
    var A = new Float64Array(p * p), B = new Float64Array(p);
    var bySchool = []; S.forEach(function (s, i) { (bySchool[s.school] = bySchool[s.school] || []).push(i); });
    bySchool.forEach(function (idx) {
      var n = idx.length, c = tau / (omega + n * tau), sx = new Float64Array(p), sy = 0;
      idx.forEach(function (i) { sy += bk[i]; for (a = 0; a < p; a++) sx[a] += S[i].x[a]; });
      idx.forEach(function (i) { var x = S[i].x; for (a = 0; a < p; a++) { B[a] += x[a] * bk[i] / omega; for (b = 0; b < p; b++) A[a * p + b] += x[a] * x[b] / omega; } });
      for (a = 0; a < p; a++) { B[a] -= c * sx[a] * sy / omega; for (b = 0; b < p; b++) A[a * p + b] -= c * sx[a] * sx[b] / omega; }
    });
    return { A: A, B: B };   // likelihood precision and precision-weighted estimate
  }

  // Posterior mean and SD for a normal likelihood (precision A, A*beta_hat = B) and an independent normal prior.
  function posterior(L, priorMean, priorSD) {
    var p = 4, A = L.A.slice(), B = L.B.slice();
    for (var a = 0; a < p; a++) { var w = 1 / (priorSD[a] * priorSD[a]); A[a * p + a] += w; B[a] += w * priorMean[a]; }
    var V = inverse(A, p); if (!V) return null;
    var m = matvec(V, B, p);
    return { mean: m, sd: Array.prototype.map.call(m, function (_, a) { return Math.sqrt(V[a * p + a]); }) };
  }

  // q: { schools, perSchool, drift (prior SDs), inflate (prior variance multiplier) }
  function replicate(q, seed, opt) {
    opt = opt || {};
    var r = rng(seed), d = simulate(q, r), sdMul = Math.sqrt(q.inflate);
    var est = { ml: [], flat: [], inf: [] }, prior = [];
    for (var k = 0; k < 3; k++) {
      var L = glsEquation(d, k); if (!L) return null;
      var V = inverse(L.A, 4), ml = matvec(V, L.B, 4);
      est.ml.push({ mean: ml, sd: [0, 1, 2, 3].map(function (a) { return Math.sqrt(V[a * 4 + a]); }) });
      est.flat.push(posterior(L, [0, 0, 0, 0], [1e4, 1e4, 1e4, 1e4]));
      var pm = TRUTH[k].map(function (v, c) { return v + (q.drift + r.normal()) * PRIOR_SD[k][c]; });
      var ps = PRIOR_SD[k].map(function (s) { return s * sdMul; });
      est.inf.push(posterior(L, pm, ps));
      prior.push({ mean: pm, sd: ps });
    }
    var out = { params: {} };
    PARAMS.forEach(function (P) {
      var row = { truth: TRUTH[P.k][P.c], prior: prior[P.k].mean[P.c], priorCenter: TRUTH[P.k][P.c] + q.drift * PRIOR_SD[P.k][P.c] };
      ["ml", "flat", "inf"].forEach(function (m) {
        var e = est[m][P.k], mu = e.mean[P.c], s = e.sd[P.c];
        row[m] = { est: mu, lo: mu - 1.96 * s, hi: mu + 1.96 * s };
      });
      out.params[P.key] = row;
    });
    if (opt.keepData) {
      // Mean curves for low (-1 SD) and high (+1 SD) SES, averaged over gender and school type.
      var nPub = Math.min(d.J - 1, Math.max(1, Math.round(PUBLIC_SHARE * d.J))), pubShare = nPub / d.J;
      var curve = function (coefs, ses) {
        return TIMES.map(function (t) {
          var v = 0;
          for (var k = 0; k < 3; k++) {
            var c = coefs[k], pi = c[0] + 0.5 * c[1] + ses * c[2] + pubShare * c[3];
            v += pi * Math.pow(t, k);
          }
          return v;
        });
      };
      var take = function (m) { return [0, 1, 2].map(function (k) { return Array.prototype.slice.call(est[m][k].mean); }); };
      out.curves = {};
      [["truth", TRUTH], ["ml", take("ml")], ["flat", take("flat")], ["inf", take("inf")]].forEach(function (pair) {
        out.curves[pair[0]] = { low: curve(pair[1], -1), high: curve(pair[1], 1) };
      });
      out.sample = d.students.filter(function (_, i) { return i % Math.max(1, Math.floor(d.students.length / 60)) === 0; })
        .map(function (s) { return { y: s.y, ses: s.x[2] }; });
      out.n = d.students.length;
    }
    return out;
  }

  return { replicate: replicate, PARAMS: PARAMS, TIMES: TIMES, TRUTH: TRUTH, PRIOR_SD: PRIOR_SD };
})();

if (typeof module !== "undefined") module.exports = GrowthSim;
