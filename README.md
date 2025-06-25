# 🌐 beachball-js

A full JavaScript/Node.js port of the `beachball.py` module from ObsPy to render earthquake focal mechanism diagrams (beachballs) on HTML canvas.

---

## 📌 Purpose

This project aims to **faithfully convert** the original Python module `beachball.py` into a modular, testable, and extensible JavaScript codebase, with support for rendering beachball diagrams using HTML5 `<canvas>`.

---

## 🚀 Features

- Full function-by-function translation from Python to JavaScript
- Written in pure ES modules for modern Node.js compatibility
- Supports rendering beachballs using canvas (`Path2D`, `context2D`)
- Built-in test framework using [Jest](https://jestjs.io/)
- Output verification against the original Python code
- Modular structure for future extensions (e.g., WebComponent support)

---

## 📁 Project Structure

```
beachball-js/
├── src/                # Main logic and helper modules
│   ├── beachball.js
│   ├── math.js
│   ├── geometry.js
│   ├── classes.js
│   └── constants.js
├── test/               # Unit tests and sample input/output
│   ├── *.test.js
│   └── data/
│       ├── input.json
│       └── expected_output.json
├── scripts/            # Python-based validation and comparison tools
├── dist/               # Minified and bundle builds for browser
├── package.json
└── README.md
```

---

## 🧪 Testing

We use **Jest** for unit testing and cross-validation with Python.

```bash
npm install
npm test
```

You can also run `scripts/compare_with_python.py` to regenerate reference outputs from the original Python code.

---

## 📚 References & Credits

This JavaScript project is a faithful adaptation of the following works:

- **Original Python Implementation:**  
  `beachball.py` by **Robert Barsch**  
  Email: `barsch@egu.eu`  
  Copyright © 2008–2012

- **Code Foundations and Algorithms Based On:**
  1. MATLAB script [`bb.m`](http://www.ceri.memphis.edu/people/olboyd/Software/Software.html)  
     by Andy Michael, Chen Ji, and Oliver Boyd  
  2. `ps_meca` from the [Generic Mapping Tools (GMT)](https://www.generic-mapping-tools.org)

- **ObsPy Python Framework:**  
  [ObsPy Development Team](https://github.com/obspy/obspy)  
  License: GNU Lesser General Public License v3  
  https://www.gnu.org/copyleft/lesser.html

This project does **not redistribute** any Python source files. All logic is **newly re-implemented** in JavaScript, while preserving mathematical and logical equivalence to the original Python source for educational and interoperability purposes.

---

## 🛠 License

This JavaScript port is released under the MIT License.  
Refer to original Python and MATLAB/GMT code for their respective licenses.

---

## 🙌 Contributions

Pull requests, bug reports, and suggestions are welcome!  
Feel free to fork this project and help make it better.

---

## 📷 Preview (Coming Soon)

Sample canvas rendering of a beachball focal mechanism will be added here.
