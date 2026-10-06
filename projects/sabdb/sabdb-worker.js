// Runs replications off the main thread and streams each one back.
importScripts("sabdb-sim.js");
onmessage = function (e) {
  var p = e.data.params, reps = e.data.reps, seed = e.data.seed;
  for (var i = 0; i < reps; i++) {
    var out = SABDBSim.replicate(p, seed + i * 7919);
    postMessage({ id: e.data.id, rep: i, reps: reps, rows: out.rows, tau: out.tau, hbTauLevel: out.hbTauLevel });
  }
};
