// Small SVG chart helpers shared by the interactive research demos.
var DemoCharts = (function () {
  var NS = "http://www.w3.org/2000/svg";
  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function txt(parent, x, y, s, attrs) {
    var t = el("text", Object.assign({ x: x, y: y }, attrs || {}), parent);
    t.textContent = s; return t;
  }
  function fmt(x, d) { return isFinite(x) ? x.toFixed(d === undefined ? 1 : d) : "–"; }
  function pct(x) { return isFinite(x) ? Math.round(x * 100) + "%" : "–"; }

  function legend(id, items) {
    var box = document.getElementById(id); box.textContent = "";
    items.forEach(function (it) {
      var s = document.createElement("span"), i = document.createElement("i");
      if (it.cls) i.className = it.cls;
      i.style.setProperty("--c", it.color); s.appendChild(i); s.appendChild(document.createTextNode(it.name)); box.appendChild(s);
    });
  }
  function niceMax(v) { var p = Math.pow(10, Math.floor(Math.log10(v))), f = v / p; return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * p; }
  function ticks(max, n) { var step = niceMax(max / n), out = []; for (var t = 0; t <= max + 1e-9; t += step) out.push(+t.toFixed(6)); return out; }
  function roundTop(x, y, w, h, r) {
    r = Math.min(r, w / 2, h);
    return "M" + x + "," + (y + h) + "V" + (y + r) + "Q" + x + "," + y + " " + (x + r) + "," + y + "H" + (x + w - r) + "Q" + (x + w) + "," + y + " " + (x + w) + "," + (y + r) + "V" + (y + h) + "Z";
  }
  function roundRight(x, y, w, h, r) {
    r = Math.min(r, h / 2, w);
    return "M" + x + "," + y + "H" + (x + w - r) + "Q" + (x + w) + "," + y + " " + (x + w) + "," + (y + r) + "V" + (y + h - r) + "Q" + (x + w) + "," + (y + h) + " " + (x + w - r) + "," + (y + h) + "H" + x + "Z";
  }
  function tipRows(tip, title, rows) {
    tip.textContent = "";
    var t = document.createElement("div"); t.className = "t-title"; t.textContent = title; tip.appendChild(t);
    rows.forEach(function (r) {
      var row = document.createElement("div"); row.className = "t-row";
      var left = document.createElement("span");
      if (r.color) { var i = document.createElement("i"); i.style.setProperty("--c", r.color); left.appendChild(i); left.appendChild(document.createTextNode(" ")); }
      left.appendChild(document.createTextNode(r.label));
      var b = document.createElement("b"); b.textContent = r.value;
      row.appendChild(left); row.appendChild(b); tip.appendChild(row);
    });
  }
  function showTip(chart, evt, title, rows) {
    var tip = chart.querySelector(".tip"); tipRows(tip, title, rows); tip.hidden = false;
    var box = chart.getBoundingClientRect(), x = evt.clientX - box.left + 14, y = evt.clientY - box.top + 14;
    if (x + tip.offsetWidth > box.width) x = Math.max(0, evt.clientX - box.left - tip.offsetWidth - 14);
    tip.style.left = x + "px"; tip.style.top = y + "px";
  }
  function hideTip(chart) { chart.querySelector(".tip").hidden = true; }
  function attachTip(node, chart, fn) {
    node.setAttribute("tabindex", "0");
    node.addEventListener("pointermove", function (e) { var t = fn(); showTip(chart, e, t[0], t[1]); });
    node.addEventListener("pointerleave", function () { hideTip(chart); });
    node.addEventListener("focus", function () {
      var r = node.getBoundingClientRect(); var t = fn();
      showTip(chart, { clientX: r.left + r.width / 2, clientY: r.top }, t[0], t[1]);
    });
    node.addEventListener("blur", function () { hideTip(chart); });
  }
  function clearSvg(chart) { var old = chart.querySelector("svg"); if (old) old.remove(); }

  return { el: el, txt: txt, fmt: fmt, pct: pct, legend: legend, niceMax: niceMax, ticks: ticks, roundTop: roundTop, roundRight: roundRight, tipRows: tipRows, showTip: showTip, hideTip: hideTip, attachTip: attachTip, clearSvg: clearSvg };
})();
