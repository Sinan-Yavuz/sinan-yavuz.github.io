// Runs replications off the main thread and streams each one back.
importScripts("mibma-sim.js");
onmessage = function (e) {
  var q = e.data.params, reps = e.data.reps, seed = e.data.seed;
  for (var i = 0; i < reps; i++) {
    var out = MIBMASim.replicate(q, seed + i * 7919, { keepData: i === 0 });
    out.id = e.data.id; out.rep = i;
    postMessage(out);
  }
};
