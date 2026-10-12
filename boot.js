// Resolves the newest published release and starts the app from the content-delivery network.
(function () {
  var root = 'https://cdn.jsdelivr.net/gh/samuelakosaonyejekwe/designermaterials', fallback = 'v1.2.0', done = false;
  function start(tag) {
    if (done) return; done = true; var base = root + '@' + tag + '/';
    var b = document.createElement('base'); b.href = base; document.head.appendChild(b);
    var l = document.createElement('link'); l.rel = 'stylesheet'; l.href = base + 'css/app.css'; document.head.appendChild(l);
    var s = document.createElement('script'); s.type = 'module'; s.src = base + 'js/app.js'; document.head.appendChild(s);
  }
  var timer = setTimeout(function () { start(fallback); }, 4000);
  fetch(root + '@main/assets.json', { cache: 'no-store' }).then(function (r) { return r.json(); }).then(function (a) { clearTimeout(timer); start(/^\d+\.\d+\.\d+$/.test(a.version) ? 'v' + a.version : fallback); }).catch(function () { clearTimeout(timer); start(fallback); });
})();
