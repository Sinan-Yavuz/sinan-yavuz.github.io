// miBMA vs mice "norm" vs complete-case analysis: in-browser simulation.
// Mirrors github.com/Sinan-Yavuz/miBMA:
//   my.bicreg: every subset of predictors is scored by BIC = n log(1 - R^2) + k log(n),
//     models outside Occam's window (BIC - min > 2 log 20) are dropped, and coefficients
//     are averaged with weights exp(-BIC / 2).
//   .norm.draw.miBMA: sigma* = sqrt(RSS / chi2_df), beta* = c + chol((X'X)^-1)' z sigma*,
//     centered on the BMA coefficients c. mice "norm" is the same draw centered on OLS.
// One incomplete variable (Y), so a single pass replaces mice's FCS iterations.

var MIBMASim = (function () {
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
    u.chisq = function (df) { return 2 * u.gamma(df / 2); };
    return u;
  }

  // Lower Cholesky factor of a p x p row-major matrix; returns null if not positive definite.
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
  function cholSolve(L, b, p) {
    var z = new Float64Array(p), x = new Float64Array(p), i, k;
    for (i = 0; i < p; i++) { var t = b[i]; for (k = 0; k < i; k++) t -= L[i * p + k] * z[k]; z[i] = t / L[i * p + i]; }
    for (i = p - 1; i >= 0; i--) { var u = z[i]; for (k = i + 1; k < p; k++) u -= L[k * p + i] * x[k]; x[i] = u / L[i * p + i]; }
    return x;
  }
  function inverse(A, p) {
    var L = chol(A, p), inv = new Float64Array(p * p);
    for (var j = 0; j < p; j++) {
      var e = new Float64Array(p); e[j] = 1;
      var col = cholSolve(L, e, p);
      for (var i = 0; i < p; i++) inv[i * p + j] = col[i];
    }
    return inv;
  }

  // Equicorrelated predictors; Y depends on the first k with equal weights scaled to the target R^2.
  // MAR: Y is missing for the cases with the largest X1 (as in the package's Example.R).
  // MCAR: Y is missing for a random subset.
  function simulate(q, r) {
    var n = q.n, p = q.p, k = Math.min(q.k, p), rho = q.rho;
    var varSignal = k + k * (k - 1) * rho;
    var b = Math.sqrt(q.r2 / (1 - q.r2) / varSignal);
    var X = [], Y = [], sr = Math.sqrt(rho), se = Math.sqrt(1 - rho);
    for (var i = 0; i < n; i++) {
      var f = r.normal(), x = new Float64Array(p), y = 0;
      for (var j = 0; j < p; j++) { x[j] = sr * f + se * r.normal(); if (j < k) y += b * x[j]; }
      X.push(x); Y.push(y + r.normal());
    }
    var order = X.map(function (_, i) { return i; });
    if (q.mech === "mcar") order.forEach(function (_, i) { var j = i + Math.floor(r() * (n - i)), t = order[i]; order[i] = order[j]; order[j] = t; });
    else order.sort(function (a, c) { return X[c][0] - X[a][0]; });
    var miss = new Uint8Array(n), nMiss = Math.round(n * q.miss);
    for (i = 0; i < nMiss; i++) miss[order[i]] = 1;
    return { X: X, Y: Y, miss: miss, p: p, k: k };
  }

  // Fit on observed rows. Returns centered quantities, OLS, and the BMA average.
  function fit(d) {
    var p = d.p, obs = [], i, j, l;
    for (i = 0; i < d.Y.length; i++) if (!d.miss[i]) obs.push(i);
    var n = obs.length, mx = new Float64Array(p), my = 0;
    obs.forEach(function (o) { my += d.Y[o]; for (j = 0; j < p; j++) mx[j] += d.X[o][j]; });
    my /= n; for (j = 0; j < p; j++) mx[j] /= n;
    var XtX = new Float64Array(p * p), Xty = new Float64Array(p), yy = 0;
    obs.forEach(function (o) {
      var x = d.X[o], yc = d.Y[o] - my;
      yy += yc * yc;
      for (j = 0; j < p; j++) {
        var xj = x[j] - mx[j]; Xty[j] += xj * yc;
        for (l = 0; l <= j; l++) XtX[j * p + l] += xj * (x[l] - mx[l]);
      }
    });
    for (j = 0; j < p; j++) for (l = j + 1; l < p; l++) XtX[j * p + l] = XtX[l * p + j];

    // Exhaustive subset search (the package uses regsubsets with nbest = 250).
    var models = [];
    for (var mask = 0; mask < (1 << p); mask++) {
      var idx = []; for (j = 0; j < p; j++) if (mask & (1 << j)) idx.push(j);
      var s = idx.length, coef = new Float64Array(p), r2 = 0;
      if (s) {
        var A = new Float64Array(s * s), bb = new Float64Array(s);
        for (var a = 0; a < s; a++) { bb[a] = Xty[idx[a]]; for (var c = 0; c < s; c++) A[a * s + c] = XtX[idx[a] * p + idx[c]]; }
        var L = chol(A, s); if (!L) continue;
        var beta = cholSolve(L, bb, s), ss = 0;
        for (a = 0; a < s; a++) { coef[idx[a]] = beta[a]; ss += beta[a] * bb[a]; }
        r2 = Math.min(ss / yy, 0.999);
      }
      models.push({ mask: mask, size: s, coef: coef, bic: n * Math.log(1 - r2) + s * Math.log(n) });
    }
    var minB = Infinity; models.forEach(function (m) { minB = Math.min(minB, m.bic); });
    var win = models.filter(function (m) { return m.bic - minB < 2 * Math.log(20); });
    var maxB = -Infinity; win.forEach(function (m) { maxB = Math.max(maxB, m.bic); });
    var wsum = 0; win.forEach(function (m) { m.w = Math.exp(-0.5 * (m.bic - maxB)); wsum += m.w; });
    var bma = new Float64Array(p), incl = new Float64Array(p), full = null;
    win.forEach(function (m) {
      m.w /= wsum;
      for (j = 0; j < p; j++) { bma[j] += m.w * m.coef[j]; if (m.mask & (1 << j)) incl[j] += m.w; }
    });
    models.forEach(function (m) { if (m.mask === (1 << p) - 1) full = m.coef; });
    return { n: n, p: p, mx: mx, my: my, XtX: XtX, XtXinv: inverse(XtX, p), yy: yy, Xty: Xty,
      ols: full, bma: bma, incl: incl, nModels: win.length };
  }

  // One imputation draw, as in .norm.draw.miBMA / mice's .norm.draw. Centered design, so
  // (X'X)^-1 is block-diagonal: 1/n for the intercept and (Xc'Xc)^-1 for the slopes.
  function draw(F, c, r) {
    var p = F.p, j, l, rss = F.yy;
    for (j = 0; j < p; j++) {
      rss -= 2 * c[j] * F.Xty[j];
      for (l = 0; l < p; l++) rss += c[j] * c[l] * F.XtX[j * p + l];
    }
    var df = Math.max(F.n - (p + 1), 1);
    var sigma = Math.sqrt(rss / r.chisq(df));
    var Lv = chol(F.XtXinv, p), z = new Float64Array(p), beta = new Float64Array(p);
    for (j = 0; j < p; j++) z[j] = r.normal();
    for (j = 0; j < p; j++) { var t = 0; for (l = 0; l <= j; l++) t += Lv[j * p + l] * z[l]; beta[j] = c[j] + t * sigma; }
    return { b0: F.my + r.normal() * sigma / Math.sqrt(F.n), beta: beta, sigma: sigma };
  }

  // Approximate 97.5% t quantile (Cornish-Fisher); exact enough for interval widths here.
  function t975(df) {
    var z = 1.959964, z3 = z * z * z, z5 = z3 * z * z;
    if (!isFinite(df) || df > 1e6) return z;
    df = Math.max(df, 2);
    return z + (z3 + z) / (4 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96 * df * df);
  }

  function replicate(q, seed, opt) {
    opt = opt || {};
    var m = opt.m || 20, r = rng(seed), d = simulate(q, r), F = fit(d), N = d.Y.length;
    var trueMean = 0;   // population mean of Y (all predictors have mean 0)
    var missIdx = []; for (var i = 0; i < N; i++) if (d.miss[i]) missIdx.push(i);

    function mi(center) {
      var Q = [], U = [], se = 0, cnt = 0, example = null;
      for (var t = 0; t < m; t++) {
        var dr = draw(F, center, r), sum = 0, ss = 0, filled = d.Y.slice();
        missIdx.forEach(function (i) {
          var x = d.X[i], yhat = dr.b0;
          for (var j = 0; j < F.p; j++) yhat += dr.beta[j] * (x[j] - F.mx[j]);
          var yimp = yhat + r.normal() * dr.sigma;
          filled[i] = yimp; se += (yimp - d.Y[i]) * (yimp - d.Y[i]); cnt++;
        });
        filled.forEach(function (y) { sum += y; }); var mean = sum / N;
        filled.forEach(function (y) { ss += (y - mean) * (y - mean); });
        Q.push(mean); U.push(ss / (N - 1) / N);
        if (t === 0) example = missIdx.map(function (i) { return filled[i]; });
      }
      var qbar = Q.reduce(function (a, b) { return a + b; }, 0) / m;
      var w = U.reduce(function (a, b) { return a + b; }, 0) / m;
      var b = Q.reduce(function (a, x) { return a + (x - qbar) * (x - qbar); }, 0) / (m - 1);
      var T = w + (1 + 1 / m) * b, rr = (1 + 1 / m) * b / w;
      var df = (m - 1) * (1 + 1 / rr) * (1 + 1 / rr), half = t975(df) * Math.sqrt(T);
      return { est: qbar, lo: qbar - half, hi: qbar + half, impRMSE: Math.sqrt(se / cnt), example: example };
    }

    var obsY = d.Y.filter(function (_, i) { return !d.miss[i]; }), no = obsY.length;
    var cm = obsY.reduce(function (a, b) { return a + b; }, 0) / no;
    var cv = obsY.reduce(function (a, y) { return a + (y - cm) * (y - cm); }, 0) / (no - 1);
    var ch = t975(no - 1) * Math.sqrt(cv / no);

    var out = {
      truth: trueMean,
      cc: { est: cm, lo: cm - ch, hi: cm + ch },
      norm: mi(F.ols),
      mibma: mi(F.bma),
      incl: Array.prototype.slice.call(F.incl),
      nModels: F.nModels, k: d.k
    };
    if (opt.keepData) {
      out.data = d.X.map(function (x, i) { return { x1: x[0], y: d.Y[i], miss: d.miss[i] }; });
      out.missIdx = missIdx;
    }
    return out;
  }

  return { replicate: replicate };
})();

if (typeof module !== "undefined") module.exports = MIBMASim;
