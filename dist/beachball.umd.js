(function (factory) {
    typeof define === 'function' && define.amd ? define(factory) :
    factory();
})((function () { 'use strict';

    //src/classes.js:
    /* ──────────────────────────────────────────────────────────
     * PrincipalAxis
     * --------------------------------------------------------*/
    class PrincipalAxis {
        /**
         * @param {number} val
         * @param {number} strike  degrees
         * @param {number} dip     degrees
         */
        constructor(val = 0, strike = 0, dip = 0) {
            this.val = val;
            this.strike = strike;
            this.dip = dip;
        }
    }

    /* ──────────────────────────────────────────────────────────
     * NodalPlane
     * --------------------------------------------------------*/
    class NodalPlane {
        /**
         * @param {number} strike degrees
         * @param {number} dip    degrees
         * @param {number} rake   degrees
         */
        constructor(strike = 0, dip = 0, rake = 0) {
            this.strike = strike;
            this.dip = dip;
            this.rake = rake;
        }
    }

    /* ──────────────────────────────────────────────────────────
     * MomentTensor
     * --------------------------------------------------------*/
    class MomentTensor {
        /**
         * Flexible constructor:
         * - (array6, expo)
         * - (matrix3x3, expo)
         * - (xx, yy, zz, xy, xz, yz, expo)  (same order as Python)
         */
        constructor(...args) {
            if (args.length === 2) {
                const a = args[0];
                this.expo = args[1];

                if (Array.isArray(a) && a.length === 6) {
                    this.mt = [
                        [a[0], a[3], a[4]],
                        [a[3], a[1], a[5]],
                        [a[4], a[5], a[2]],
                    ];
                } else if (Array.isArray(a) && a.length === 3 && Array.isArray(a[0])) {
                    // Assume 3×3 full matrix
                    if (a.length !== 3 || a[0].length !== 3) {
                        throw new TypeError('Full matrix must be 3×3.');
                    }
                    this.mt = a;
                } else {
                    throw new TypeError('Invalid parameter size for MomentTensor.');
                }
            } else if (args.length === 7) {
                const [xx, yy, zz, xy, xz, yz, expo] = args;
                this.mt = [
                    [xx, xy, xz],
                    [xy, yy, yz],
                    [xz, yz, zz],
                ];
                this.expo = expo;
            } else {
                throw new TypeError('Wrong number of parameters for MomentTensor.');
            }
        }

        /* ---------- derived helpers ---------- */

        get mtNormalized() {
            const norm = this.mtFrobeniusNorm();
            return this.mt.map(row => row.map(v => v / norm));
        }

        mtFrobeniusNorm() {
            return Math.sqrt(
                this.mt.reduce((acc, row) => acc + row.reduce((a, v) => a + v * v, 0), 0)
            );
        }

        /** Returns a new normalized MomentTensor */
        get normalized() {
            return new MomentTensor(this.mtNormalized, this.expo);
        }

        /* ---------- component getters ---------- */
        get xx() { return this.mt[0][0]; }
        get xy() { return this.mt[0][1]; }
        get xz() { return this.mt[0][2]; }
        get yy() { return this.mt[1][1]; }
        get yz() { return this.mt[1][2]; }
        get zz() { return this.mt[2][2]; }
    }

    /**
     * eigenSym3(A, order = 'asc')
     * @param {number[][]} A    3×3 symmetric matrix
     * @param {'asc' | 'desc'} order  Sort order for eigenvalues and eigenvectors
     * @returns {{values: number[], vectors: number[][]}} Sorted eigenvalues and corresponding eigenvectors
     */
    function eigenSym3(A, order = 'asc') {
        // A is a symmetric 3×3 matrix: [[a11,a12,a13],[a12,a22,a23],[a13,a23,a33]]

        /**
         * Normalize vector sign so the element of largest absolute value is non-negative
         */
        function normalizeSign(v) {
            // find index of max absolute component
            const maxIdx = v.reduce((iMax, x, i, arr) =>
                Math.abs(x) > Math.abs(arr[iMax]) ? i : iMax, 0);
            return v[maxIdx] < 0 ? v.map(x => -x) : v;
        }

        /**
         * Sort eigenvalues and eigenvectors according to the specified order
         */
        function sort(values, vectors) {
            // For ascending: [λ_min, λ_mid, λ_max] => idx = [2,1,0]
            // For descending: [λ_max, λ_mid, λ_min] => idx = [0,1,2]
            const idx = order === 'asc' ? [2, 1, 0] : [0, 1, 2];
            const sortedValues = idx.map(i => values[i]);
            const sortedVectors = idx.map(i => normalizeSign(vectors[i]));
            return { values: sortedValues, vectors: sortedVectors };
        }

        // 1) Scale matrix to avoid overflow
        let s = 0;
        for (let i = 0; i < 3; ++i) {
            for (let j = i; j < 3; ++j) {
                s = Math.max(s, Math.abs(A[i][j]));
            }
        }
        if (s === 0) {
            // Zero matrix -> zero eigenvalues and standard basis eigenvectors
            return sort([0, 0, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        }
        const M = A.map(row => row.map(val => val / s));

        // 2) Compute eigenvalues using the analytical method (Kopp)
        const [m11, m12, m13] = [M[0][0], M[0][1], M[0][2]];
        const [m22, m23, m33] = [M[1][1], M[1][2], M[2][2]];
        const p1 = m12 * m12 + m13 * m13 + m23 * m23;

        if (p1 === 0) {
            // Diagonal matrix case
            const diagVals = [m11, m22, m33].map(v => v * s);
            return sort(diagVals, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        }
        const q = (m11 + m22 + m33) / 3;
        const p2 = (m11 - q) ** 2 + (m22 - q) ** 2 + (m33 - q) ** 2 + 2 * p1;
        const p = Math.sqrt(p2 / 6);

        // Compute B = (M - q I) / p
        const B11 = (m11 - q) / p, B12 = m12 / p, B13 = m13 / p;
        const B22 = (m22 - q) / p, B23 = m23 / p, B33 = (m33 - q) / p;

        // Compute the determinant-like quantity r
        const r = (
            B11 * B22 * B33 + 2 * B12 * B13 * B23
            - B11 * B23 * B23 - B22 * B13 * B13 - B33 * B12 * B12
        ) / 2;
        const rClamped = Math.max(-1, Math.min(1, r));
        const phi = Math.acos(rClamped) / 3;

        // Eigenvalues of the normalized matrix
        const eig1 = q + 2 * p * Math.cos(phi);               // largest
        const eig3 = q + 2 * p * Math.cos(phi + 2 * Math.PI / 3); // smallest
        const eig2 = 3 * q - eig1 - eig3;                     // mid
        const values = [eig1, eig2, eig3].map(v => v * s);

        // 3) Compute eigenvectors via cross-product method
        function eigenvectorFor(lambda) {
            const m = [
                [M[0][0] - lambda, M[0][1], M[0][2]],
                [M[1][0], M[1][1] - lambda, M[1][2]],
                [M[2][0], M[2][1], M[2][2] - lambda]
            ];
            let v = [
                m[0][1] * m[1][2] - m[0][2] * m[1][1],
                m[0][2] * m[1][0] - m[0][0] * m[1][2],
                m[0][0] * m[1][1] - m[0][1] * m[1][0]
            ];
            let norm = Math.hypot(...v);
            if (norm === 0) {
                v = [
                    m[0][1] * m[2][2] - m[0][2] * m[2][1],
                    m[0][2] * m[2][0] - m[0][0] * m[2][2],
                    m[0][0] * m[2][1] - m[0][1] * m[2][0]
                ];
                norm = Math.hypot(...v);
                if (norm === 0) { v = [0, 0, 0]; norm = 1; }
            }
            return v.map(x => x / norm);
        }

        const rawVectors = values.map(val => eigenvectorFor(val / s));

        // 4) Return sorted (and LAPACK‐style normalized) eigenpairs
        return sort(values, rawVectors);
    }

    // geometry.js:

    /**
     * Finds strike and dip of a plane given its normal-vector components (n, e, u).
     * Mirrors the exact logic of the original Python implementation.
     *
     * @param {number} n - North component of the normal vector
     * @param {number} e - East  component of the normal vector
     * @param {number} u - Up    component of the normal vector
     * @returns {[number, number]} [strikeDeg, dipDeg] in degrees
     */
    function strikeDip(n, e, u) {
        const r2d = 180 / Math.PI;

        // Ensure upward-pointing normal
        if (u < 0) {
            if (n != 0) n = -n;
            if (e != 0) e = -e;
            u = -u;
        }

        // Strike computation
        let strike = Math.atan2(e, n) * r2d;
        strike -= 90;

        while (strike >= 360) strike -= 360;
        while (strike < 0) strike += 360;

        // Dip computation
        const x = Math.sqrt(n * n + e * e);
        const dip = Math.atan2(x, u) * r2d;

        return [strike, dip];
    }

    /**
     * Generates a closed polygon path from x/y point arrays with optional scaling
     * (res) and translation (xy). The returned structure is framework-agnostic:
     *   { vertices: [[x0, y0], [x1, y1], ...],
     *     codes:    ['M', 'L', 'L', ..., 'Z'] }
     *
     * Canvas usage example in browser:
     *   const p = xy2patch([...], [...], [sx, sy], [cx, cy]);
     *   const path = new Path2D();
     *   p.vertices.forEach(([vx, vy], i) => {
     *     if (p.codes[i] === 'M') path.moveTo(vx, vy);
     *     else if (p.codes[i] === 'L') path.lineTo(vx, vy);
     *     else if (p.codes[i] === 'Z') path.closePath();
     *   });
     *   ctx.fill(path);
     *
     * @param {number[]} x   – array of x-coordinates
     * @param {number[]} y   – array of y-coordinates
     * @param {number|number[]} res – scale factor(s); scalar or [sx, sy]
     * @param {number[]} xy  – translation [cx, cy]
     * @returns {{vertices: number[][], codes: string[]}}
     */
    function xy2patch(x, y, res, xy) {
        // ── 1. Normalise `res` to a two-component array ──
        try {
            // If 'res' is a single primitive (number), it has no 'length' property
            if (res.length === undefined) {
                throw new TypeError();          // Mirrors Python’s TypeError in the assert
            }

            // If iterable, it must contain exactly two elements
            if (res.length !== 2) {
                throw new Error('res must contain exactly two elements'); // Like AssertionError
            }
        } catch (err) {
            // Only handle the TypeError that signals a single-value resolution
            if (err instanceof TypeError) {
                res = [res, res];               // Same resolution in both axes → circle
            } else {
                throw err;                      // Re-throw any other errors
            }
        }

        // ── 2. Transform every (x, y) into the patch coordinate system ──
        const vertices = x.map((xi, i) => [
            xi * res[0] + xy[0],            // Scale & translate X
            y[i] * res[1] + xy[1]           // Scale & translate Y
        ]);

        // ── 3. Return **data only**, no drawing primitives ──
        return {
            vertices,                       // [[x0,y0], [x1,y1], …] – closed polygon
            res,                            // Saved for potential re-scaling
            center: xy                      // Translation vector, useful for later
        };                                 // equivalent to PathPatch
    }

    /**
     * Calculates strike, dip, and rake of the auxiliary (second) plane given the
     * first plane's strike (s1), dip (d1), and rake (r1). All angles in degrees.
     *
     * @param {number} s1  Strike of first plane (deg)
     * @param {number} d1  Dip    of first plane (deg)
     * @param {number} r1  Rake   of first plane (deg)
     * @returns {[number, number, number]} [strike, dip, rake] of second plane (deg)
     */
    function auxPlane(s1, d1, r1) {
        const r2d = 180 / Math.PI;

        const z = (s1 + 90) / r2d;
        const z2 = d1 / r2d;
        const z3 = r1 / r2d;

        // Slick vector in plane 1
        const sl1 = -Math.cos(z3) * Math.cos(z) - Math.sin(z3) * Math.sin(z) * Math.cos(z2);
        const sl2 = Math.cos(z3) * Math.sin(z) - Math.sin(z3) * Math.cos(z) * Math.cos(z2);
        const sl3 = Math.sin(z3) * Math.sin(z2);

        const [strike, dip] = strikeDip(sl2, sl1, sl3);

        // Normal vector to plane 1
        const n1 = Math.sin(z) * Math.sin(z2);
        const n2 = Math.cos(z) * Math.sin(z2);

        // Strike vector of plane 2
        const h1 = -sl2;
        const h2 = sl1;

        // Rake calculation
        let zz = (h1 * n1 + h2 * n2) / Math.sqrt(h1 * h1 + h2 * h2);

        const float64epsilon = 2.2204460492503131e-16;
        if (Math.abs(zz) > 1 && Math.abs(zz) < 1 + 100 * float64epsilon) {
            zz = Math.sign(zz);
        }

        const zAngle = Math.acos(zz);
        const rake = sl3 > 0 ? zAngle * r2d : -zAngle * r2d;

        return [strike, dip, rake];
    }

    /**
     * Calculates three angles (ft, fd, fl) from vectors an and bn.
     * Used internally by mt2plane. All angles in degrees.
     *
     * @param {number[]} an - vector a (length 3)
     * @param {number[]} bn - vector b (length 3)
     * @returns {[number, number, number]} [ft, fd, fl] in degrees
     */
    function tdl(an, bn) {
        let [xn, yn, zn] = an;
        const [xe, ye, ze] = bn;

        const aaa = 1.0 / 1_000_000;
        const con = 57.2957795; // rad → deg conversion factor

        let fd, ft, fl;

        if (Math.abs(zn) < aaa) {
            fd = 90.0;

            let axn = Math.abs(xn);
            if (axn > 1.0) axn = 1.0;
            ft = Math.asin(axn) * con;

            const st = -xn;
            const ct = yn;

            if (st >= 0 && ct < 0) ft = 180.0 - ft;
            if (st < 0 && ct <= 0) ft = 180.0 + ft;
            if (st < 0 && ct > 0) ft = 360.0 - ft;

            fl = Math.asin(Math.abs(ze)) * con;
            const sl = -ze;
            let cl;
            if (Math.abs(xn) < aaa) cl = xe / yn;
            else cl = -ye / xn;

            if (sl >= 0 && cl < 0) fl = 180.0 - fl;
            if (sl < 0 && cl <= 0) fl -= 180.0;
            if (sl < 0 && cl > 0) fl = -fl;

        } else {
            if (-zn > 1.0) zn = -1.0;

            const fdh = Math.acos(-zn);
            fd = fdh * con;

            const sd = Math.sin(fdh);
            if (sd === 0) return; // singular

            const st = -xn / sd;
            const ct = yn / sd;

            let sx = Math.abs(st);
            if (sx > 1.0) sx = 1.0;
            ft = Math.asin(sx) * con;

            if (st >= 0 && ct < 0) ft = 180.0 - ft;
            if (st < 0 && ct <= 0) ft = 180.0 + ft;
            if (st < 0 && ct > 0) ft = 360.0 - ft;

            let sl = -ze / sd;
            sx = Math.abs(sl);
            if (sx > 1.0) sx = 1.0;
            fl = Math.asin(sx) * con;

            let cl;
            if (st === 0) {
                cl = xe / ct;
            } else {
                const xxx = yn * zn * ze / (sd * sd) + ye;
                cl = -sd * xxx / xn;
                if (ct === 0) cl = ye / st;
            }

            if (sl >= 0 && cl < 0) fl = 180.0 - fl;
            if (sl < 0 && cl <= 0) fl -= 180.0;
            if (sl < 0 && cl > 0) fl = -fl;
        }

        return [ft, fd, fl];
    }

    // mt2plane.js :

    /**
     * Computes a NodalPlane object from a given MomentTensor.
     *
     * @param {MomentTensor} mt - Moment tensor instance
     * @returns {NodalPlane}
     */
    function mt2plane(mt) {
        // 1) eigen-decompose in ascending order (λmin→λmid→λmax)
        const { values: d0, vectors: ev } = eigenSym3(mt.mt, 'asc');

        // 2) build a temporary v0[row][col]
        const v0 = [
            [ev[0][0], ev[1][0], ev[2][0]],
            [ev[0][1], ev[1][1], ev[2][1]],
            [ev[0][2], ev[1][2], ev[2][2]],
        ];

        // 3) “Python trick” to re-order and re-sign the eigenvector matrix
        const d = [d0[1], d0[2], d0[0]];   // [mid, max, min]
        let v = [
            [v0[1][1], -v0[1][2], -v0[1][0]],  // row 0
            [v0[2][1], -v0[2][2], -v0[2][0]],  // row 1
            [-v0[0][1], v0[0][2], v0[0][0]], // row 2
        ];

        // — now force each column j to have v[0][j] ≥ 0, just like Python’s mapping does —
        for (let j = 0; j < 3; j++) {
            if (v[0][j] < 0) {
                v[0][j] *= -1;
                v[1][j] *= -1;
                v[2][j] *= -1;
            }
        }

        /* 4) ادامهٔ محاسبات (ae / an …) بدون تغییر */
        const imax = d.indexOf(Math.max(...d));
        const imin = d.indexOf(Math.min(...d));

        const ae = [
            (v[0][imax] + v[0][imin]) / Math.SQRT2,
            (v[1][imax] + v[1][imin]) / Math.SQRT2,
            (v[2][imax] + v[2][imin]) / Math.SQRT2,
        ];
        const an = [
            (v[0][imax] - v[0][imin]) / Math.SQRT2,
            (v[1][imax] - v[1][imin]) / Math.SQRT2,
            (v[2][imax] - v[2][imin]) / Math.SQRT2,
        ];

        const norm = u => Math.hypot(...u);
        const aeN = ae.map(x => x / norm(ae));
        const anN = norm(an) === 0 ? [NaN, NaN, NaN] : an.map(x => x / norm(an));

        const flip = anN[2] > 0;
        const an1 = flip ? anN.map(x => -x) : anN;
        const ae1 = flip ? aeN.map(x => -x) : aeN;

        const [ft, fd, fl] = tdl(an1, ae1);
        return new NodalPlane(360 - ft, fd, 180 - fl);
    }

    //trigHelpers.js:
    //trigonometry Helpers

    // Angle conversions
    const D2R = Math.PI / 180;  // degrees → radians
    const R2D = 180 / Math.PI;  // radians → degrees

    // Small tolerance value
    const EPSILON = 1e-5;



    /**
    * Notes:
    * - If either input is an array, the other is broadcast to the same length.
    * - Throws if array lengths mismatch.
     *
     * @param {number | number[]} th angle(s) in radians
     * @param {number | number[]} r radius/radii
     * @returns {[number | number[], number | number[]]} (same shape as inputs)
     */
    function pol2cart(th, r) {
        const isArrayTH = Array.isArray(th);
        const isArrayR = Array.isArray(r);

        // Helper for scalar calculations
        const calc = (theta, radius) => {
            const xVal = radius * Math.cos(theta);
            const yVal = radius * Math.sin(theta);
            return [xVal, yVal];
        };

        // Case 1: both scalars
        if (!isArrayTH && !isArrayR) {
            return calc(th, r);
        }

        // Case 2: broadcast arrays
        const thArr = isArrayTH ? th : Array(r.length).fill(th);
        const rArr = isArrayR ? r : Array(th.length).fill(r);

        if (thArr.length !== rArr.length) {
            throw new Error('pol2cart: th and r must have the same length when arrays');
        }

        const xs = [];
        const ys = [];
        for (let i = 0; i < thArr.length; i += 1) {
            const [xi, yi] = calc(thArr[i], rArr[i]);
            xs.push(xi);
            ys.push(yi);
        }
        return [xs, ys];
    }

    /**
     * Calculates principal axes T, N, P from a moment tensor.
     *
     * @param {MomentTensor} mt
     * @returns {[PrincipalAxis, PrincipalAxis, PrincipalAxis]}
     */
    function mt2axes(mt) {
        // 1) eigen-decompose (any order)
        const { values: d, vectors: ev } = eigenSym3(mt.mt, 'asc');  // λmin→λmax

        // 2) build v[row][col]
        const v = [
            [ev[0][0], ev[1][0], ev[2][0]],
            [ev[0][1], ev[1][1], ev[2][1]],
            [ev[0][2], ev[1][2], ev[2][2]],
        ];

        // 3) compute plunge & azimuth arrays
        let pl = [
            Math.asin(-v[0][0]),
            Math.asin(-v[0][1]),
            Math.asin(-v[0][2]),
        ];
        let az = [
            Math.atan2(v[2][0], -v[1][0]),
            Math.atan2(v[2][1], -v[1][1]),
            Math.atan2(v[2][2], -v[1][2]),
        ];

        for (let i = 0; i < 3; i++) {
            if (pl[i] <= 0) { pl[i] = -pl[i]; az[i] += Math.PI; }
            if (az[i] < 0) az[i] += 2 * Math.PI;
            if (az[i] > 2 * Math.PI) az[i] -= 2 * Math.PI;
        }

        pl = pl.map(r => r * R2D);
        az = az.map(r => r * R2D);

        // 4) find indices of min, mid, max in d
        const sortedIdx = [0, 1, 2].sort((i, j) => d[i] - d[j]);
        const [iMin, iMid, iMax] = sortedIdx;

        // 5) construct axes in correct order
        const T = new PrincipalAxis(d[iMax], az[iMax], pl[iMax]);
        const N = new PrincipalAxis(d[iMid], az[iMid], pl[iMid]);
        const P = new PrincipalAxis(d[iMin], az[iMin], pl[iMin]);

        return [T, N, P];
    }

    /**
     * Build the geometric description of a double-couple focal-mechanism
     * “beach-ball” suitable for later rendering on any 2-D backend
     * (Canvas, SVG, WebGL …).
     *
     * No drawing is performed here; the function returns *data only*.
     *
     * @param {Object}        np1              – First nodal-plane object
     * @param {number}        np1.strike       – Strike  (degrees, 0–360)
     * @param {number}        np1.dip          – Dip     (degrees, 0–90)
     * @param {number}        np1.rake         – Rake    (degrees, −180–180 or 0–360)
     * @param {number} [size = 200]            – Output square size in pixels
     * @param {[number,number]} [xy = [0,0]]   – Centre of the plot in user units
     * @param {number|[number,number]} [width = 200]
     *        If a single number, circle radius in pixels.
     *        If a pair [xw, yw], ellipse half-widths along *x* and *y*.
     *
     * @returns {[string[], Object[]]}
     *   • **colours** — two-element array `['b','w']` (pressure colour first).  
     *   • **patches** — array of objects produced by `xy2patch()`, each with:
     *     ─ `vertices` : `[x,y][]` closed polygon (screen units)  
     *     ─ `res`      : `[rx,ry]` per-axis resolution scaling (unitless)  
     *     ─ `center`   : `[cx,cy]` translation applied when creating vertices
     */
    function plotDC(np1, size = 200, xy = [0, 0], width = 200) {
        // ─── 1. Check whether 'width' indicates a Circle (one value) or an Ellipse (two values) ───
        try {
            // If 'width' is a single primitive (e.g. number), it has no 'length' property
            if (width.length === undefined) {
                throw new TypeError();          // Mirrors Python’s TypeError in the assert
            }

            // If iterable, it must contain exactly two elements
            if (width.length !== 2) {
                throw new Error('width must contain exactly two elements'); // Like AssertionError
            }
        } catch (err) {
            // Handle only the TypeError that signals a single-value width
            if (err instanceof TypeError) {
                width = [width, width];         // Same width in both directions → circle
            } else {
                throw err;                      // Re-throw any other errors
            }
        }

        // ─── 2. Normalise rake angle to the 0–180° interval ───
        let s_1 = np1.strike;       // Strike
        let d_1 = np1.dip;          // Dip
        let r_1 = np1.rake;         // Rake

        let m = 0;                  // Flag: 1 if rake was adjusted

        // If rake exceeds 180°, wrap it down into the 0–180° range
        if (r_1 > 180) {
            r_1 -= 180;
            m = 1;                    // Original rake was out of range
        }

        // If rake is negative, shift it up into the 0–180° range
        if (r_1 < 0) {
            r_1 += 180;
            m = 1;                    // Original rake was out of range
        }


        // Get azimuth and dip of the second plane
        let [s_2, d_2, _r_2] = auxPlane(s_1, d_1, r_1); 

        // ─── 3. Prepare half-size and clamp dip angles below 90° ───
        const d = size / 2;                 // Half of the canvas size (radius in pixels)

        // Dip angles of exactly 90° cause singularities; cap them at 89.9999°
        if (d_1 >= 90) d_1 = 89.9999;
        if (d_2 >= 90) d_2 = 89.9999;


        // ─── 4. Create a phi array (0 … π, step 0.01) ───
        // Note: Math.PI is not an exact multiple of 0.01, so the last value is < π
        const phi = [];
        for (let t = 0; t < Math.PI; t += 0.01) {
            phi.push(t);
        }

        // ─── 5. Compute l1 and l2 for the two nodal planes ───
        const ninetySquared = 90 * 90;
        const k1 = 90 - d_1;               // Constant term for first plane
        const k2 = 90 - d_2;               // Constant term for second plane

        const l1 = phi.map(t => {
            // numerator = (90 − d1)²
            // denominator = sin²φ + cos²φ · (90 − d1)² / 90²
            const numerator = k1 * k1;
            const denominator =
                Math.sin(t) ** 2 +
                (Math.cos(t) ** 2) * (k1 * k1) / ninetySquared;
            return Math.sqrt(numerator / denominator);
        });

        const l2 = phi.map(t => {
            const numerator = k2 * k2;
            const denominator =
                Math.sin(t) ** 2 +
                (Math.cos(t) ** 2) * (k2 * k2) / ninetySquared;
            return Math.sqrt(numerator / denominator);
        });

        // ─── 6. Build polygon paths for tension (white) and pressure (black) lobes ───
        const collect = [];

        // Iterate first over the tension side, then over the pressure side
        for (const m_ of [(m + 1) % 2, m]) {
            let inc = 1;

            // Segment on first nodal great-circle
            const [x_1, y_1] = pol2cart(
                phi.map(t => t + s_1 * D2R),   // φ shifted by strike of plane-1
                l1                              // corresponding radii
            );

            let th1, th2, xs_1, ys_1, x_2, y_2;

            if (m_ === 1) {
                // ----- Tension lobe -----
                let lo = s_1 - 180;
                let hi = s_2;
                if (lo > hi) inc = -1;

                th1 = range(s_1 - 180, s_2, inc);           // mimic np.arange
                [xs_1, ys_1] = pol2cart(
                    th1.map(t => t * D2R),
                    Array(th1.length).fill(90)
                );

                [x_2, y_2] = pol2cart(
                    phi.map(t => t + s_2 * D2R),
                    l2
                );

                th2 = range(s_2 + 180, s_1, -inc);
            } else {
                // ----- Pressure lobe -----
                let hi = s_1 - 180;
                let lo = s_2 - 180;
                if (lo > hi) inc = -1;

                th1 = range(hi, lo, -inc);
                [xs_1, ys_1] = pol2cart(
                    th1.map(t => t * D2R),
                    Array(th1.length).fill(90)
                );

                [x_2, y_2] = pol2cart(
                    phi.map(t => t + s_2 * D2R),
                    l2
                );
                x_2.reverse();                     // reverse to preserve winding
                y_2.reverse();

                th2 = range(s_2, s_1, inc);
            }

            // Closing arc on the auxiliary great-circle
            const [xs_2, ys_2] = pol2cart(
                th2.map(t => t * D2R),
                Array(th2.length).fill(90)
            );

            // Concatenate all pieces into one polygon
            const x = [...x_1, ...xs_1, ...x_2, ...xs_2];
            const y = [...y_1, ...ys_1, ...y_2, ...ys_2];

            // Scale from degrees (radius 90) to pixel units (radius d)
            const xScaled = x.map(v => v * d / 90);
            const yScaled = y.map(v => v * d / 90);

            // Pixel resolution relative to canvas size
            const res = width.map(v => v / size);

            // Convert (x,y) arrays to a drawable patch and store
            collect.push(xy2patch(yScaled, xScaled, res, xy));
        }

        // Helper: numpy-like arange
        function range(start, stop, step = 1) {
            const out = [];
            if (step > 0) {
                for (let v = start; v < stop; v += step) out.push(v);
            } else {
                for (let v = start; v > stop; v += step) out.push(v);
            }
            return out;
        }

        return [['b', 'w'], collect];   // colour order + collection of polygons
    }

    // plotMT.js:
     
    /**
     * Generate patch objects and corresponding fill colors for a moment-tensor beachball.
     *
     * @param {PrincipalAxis} T             Tension principal axis.
     * @param {PrincipalAxis} N             Neutral/principal axis.
     * @param {PrincipalAxis} P             Pressure principal axis.
     * @param {number}        size          Diameter of the beachball in pixels.
     * @param {boolean}       plotZeroTrace If true, force zero isotropic component.
     * @param {number}        x0            X-coordinate offset of the center.
     * @param {number}        y0            Y-coordinate offset of the center.
     * @param {[number,number]} xy          Center position as [x, y].
     * @param {[number,number]|number} width Ellipse radii: single value or [width, height].
     * @returns {[string[], object[]]}       A tuple of:
     *   - colors: array of fill-color strings (e.g. `"b"`, `"w"`),  
     *   - patches: array of patch-definition objects for drawing.
     */
    function plotMT(T, N, P, size = 200, plot_zerotrace = true, x0 = 0, y0 = 0, xy = [0, 0], width = 200) {

        /* ---- check if one or two widths are specified (Circle / Ellipse) */
        /* from matplotlib import patches  ← irrelevant in JS               */
        try {
            if (!Array.isArray(width) || width.length !== 2) throw new TypeError();
        } catch (e) {
            width = [width, width];
        }

        /* ---- Python lists replaced with JS arrays of same length -------- */
        const collect = [];
        const colors = [];
        const res = width.map(value => value / size);
        let b = 1;
        let big_iso = 0;
        let j = 1;
        let j2 = 0;
        let j3 = 0;
        let n = 0;

        const azi = Array.from({ length: 3 }, () => new Array(2).fill(0));
        const x = new Array(400).fill(0);
        const y = new Array(400).fill(0);
        const x2 = new Array(400).fill(0);
        const y2 = new Array(400).fill(0);
        const x3 = new Array(400).fill(0);
        const y3 = new Array(400).fill(0);
        const xp1 = new Array(800).fill(0);
        const yp1 = new Array(800).fill(0);
        const xp2 = new Array(400).fill(0);
        const yp2 = new Array(400).fill(0);

        /* ---- a, p, v arrays -------------------------------------------- */
        const a = [T.strike, N.strike, P.strike];
        const p = [T.dip, N.dip, P.dip];
        const v = [T.val, N.val, P.val];

        let vi = (v[0] + v[1] + v[2]) / 3.0;
        for (let i = 0; i < 3; i++) {
            v[i] = v[i] - vi;
        }

        const radius_size = size * 0.5;

        /* ---- pure implosion-explosion check ----------------------------- */
        if (Math.abs(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) < EPSILON) {
            if (vi > 0.0) {
                const cir = { type: "ellipse", xy, width: width[0], height: width[1] };
                collect.push(cir);
                colors.push("b");
            }
            if (vi < 0.0) {
                const cir = { type: "ellipse", xy, width: width[0], height: width[1] };
                collect.push(cir);
                colors.push("w");
            }
            return [colors, collect];
        }

        /* ---- choose d / m ---------------------------------------------- */
        let d, m;
        if (Math.abs(v[0]) >= Math.abs(v[2])) { d = 0; m = 2; }
        else { d = 2; m = 0; }

        if (plot_zerotrace) vi = 0.0;

        const f = -v[1] / v[d];
        const iso = vi / v[d];

        // # Cliff Frohlich, Seismological Research letters,
        // # Vol 7, Number 1, January-February, 1996
        // # Unless the isotropic parameter lies in the range
        // # between -1 and 1 - f there will be no nodes whatsoever
        /* ---- Frohlich node existence test ------------------------------ */
        if (iso < -1) {
            collect.push({ type: "ellipse", xy, width: width[0], height: width[1] });
            colors.push("w");
            return [colors, collect];
        } else if (iso > 1 - f) {
            collect.push({ type: "ellipse", xy, width: width[0], height: width[1] });
            colors.push("b");
            return [colors, collect];
        }

        /* ---- lots of trig shorthands ----------------------------------- */
        const spd = Math.sin(p[d] * D2R);
        const cpd = Math.cos(p[d] * D2R);
        const spb = Math.sin(p[b] * D2R);
        const cpb = Math.cos(p[b] * D2R);
        const spm = Math.sin(p[m] * D2R);
        const cpm = Math.cos(p[m] * D2R);
        const sad = Math.sin(a[d] * D2R);
        const cad = Math.cos(a[d] * D2R);
        const sab = Math.sin(a[b] * D2R);
        const cab = Math.cos(a[b] * D2R);
        const sam = Math.sin(a[m] * D2R);
        const cam = Math.cos(a[m] * D2R);

        /* ---- main 0-to-359 loop ---------------------------------------- */
        let azp = 0.0;
        for (let i = 0; i < 360; i++) {
            const fir = i * D2R;
            const s2alphan = (2.0 + 2.0 * iso) /
                (3.0 + (1.0 - 2.0 * f) * Math.cos(2.0 * fir));

            if (s2alphan > 1.0) {
                big_iso += 1;
            } else {
                const alphan = Math.asin(Math.sqrt(s2alphan));
                const sfi = Math.sin(fir);
                const cfi = Math.cos(fir);
                const san = Math.sin(alphan);
                const can = Math.cos(alphan);

                const xz = can * spd + san * sfi * spb + san * cfi * spm;
                const xn = can * cpd * cad + san * sfi * cpb * cab +
                    san * cfi * cpm * cam;
                const xe = can * cpd * sad + san * sfi * cpb * sab +
                    san * cfi * cpm * sam;

                let takeoff, az;
                if (Math.abs(xn) < EPSILON && Math.abs(xe) < EPSILON) {
                    takeoff = 0.0;
                    az = 0.0;
                } else {
                    az = Math.atan2(xe, xn);
                    if (az < 0.0) az += Math.PI * 2.0;
                    takeoff = Math.acos(
                        xz / Math.sqrt(xz * xz + xn * xn + xe * xe)
                    );
                }

                if (takeoff > Math.PI / 2.0) {
                    takeoff = Math.PI - takeoff;
                    az += Math.PI;
                    if (az > Math.PI * 2.0) az -= Math.PI * 2.0;
                }

                const r = Math.SQRT2 * Math.sin(takeoff / 2.0);
                const si = Math.sin(az);
                const co = Math.cos(az);

                if (i === 0) {
                    azi[i][0] = az;
                    x[i] = x0 + radius_size * r * si;
                    y[i] = y0 + radius_size * r * co;
                    azp = az;
                } else {
                    if (Math.abs(Math.abs(az - azp) - Math.PI) < D2R * 10.0) {
                        azi[n][1] = azp;
                        n += 1;
                        azi[n][0] = az;
                    }
                    if (Math.abs(Math.abs(az - azp) - Math.PI * 2.0) < D2R * 2.0) {
                        if (azp < az) azi[n][0] += Math.PI * 2.0;
                        else azi[n][0] -= Math.PI * 2.0;
                    }
                    if (n === 0) {
                        x[j] = x0 + radius_size * r * si;
                        y[j] = y0 + radius_size * r * co;
                        j += 1;
                    } else if (n === 1) {
                        x2[j2] = x0 + radius_size * r * si;
                        y2[j2] = y0 + radius_size * r * co;
                        j2 += 1;
                    } else if (n === 2) {
                        x3[j3] = x0 + radius_size * r * si;
                        y3[j3] = y0 + radius_size * r * co;
                        j3 += 1;
                    }
                    azp = az;
                }
            }
        }
        azi[n][1] = azp;

        /* ---- color selections based on v[1] ----------------------------- */
        let rgb1, rgb2;
        if (v[1] < 0.0) {
            rgb1 = "b";
            rgb2 = "w";
        } else {
            rgb1 = "w";
            rgb2 = "b";
        }

        /* ---- background ellipse ---------------------------------------- */
        const cir = { type: "ellipse", xy, width: width[0], height: width[1] };
        collect.push(cir);
        colors.push(rgb2);

        /* ---- n == 0 case ------------------------------------------------ */
        if (n === 0) {
            collect.push(xy2patch(x.slice(0, 360), y.slice(0, 360), res, xy)); 
            colors.push(rgb1);
            return [colors, collect];
        }

        /* ---- n == 1 case ------------------------------------------------ */
        if (n === 1) {
            /* copy x,y into xp1,yp1 */
            let i;
            for (i = 0; i < j; i++) {
                xp1[i] = x[i];
                yp1[i] = y[i];
            }
            i -= 1;

            /* unwrap azi[0][0] w.r.t. azi[0][1] */
            if (azi[0][0] - azi[0][1] > Math.PI) {
                azi[0][0] -= Math.PI * 2.0;
            } else if (azi[0][1] - azi[0][0] > Math.PI) {
                azi[0][0] += Math.PI * 2.0;
            }

            /* add arc points between azi[0][1] and azi[0][0] */
            if (azi[0][0] < azi[0][1]) {
                let az = azi[0][1] - D2R;
                while (az > azi[0][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp1[i] = x0 + radius_size * si;
                    yp1[i] = y0 + radius_size * co;
                    i += 1;
                    az -= D2R;
                }
            } else {
                let az = azi[0][1] + D2R;
                while (az < azi[0][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp1[i] = x0 + radius_size * si;
                    yp1[i] = y0 + radius_size * co;
                    i += 1;
                    az += D2R;
                }
            }

            collect.push(xy2patch(xp1.slice(0, i), yp1.slice(0, i), res, xy));
            colors.push(rgb1);

            /* copy x2,y2 into xp2,yp2 */
            for (i = 0; i < j2; i++) {
                xp2[i] = x2[i];
                yp2[i] = y2[i];
            }
            i -= 1;

            /* unwrap azi[1][0] w.r.t azi[1][1] */
            if (azi[1][0] - azi[1][1] > Math.PI) {
                azi[1][0] -= Math.PI * 2.0;
            } else if (azi[1][1] - azi[1][0] > Math.PI) {
                azi[1][0] += Math.PI * 2.0;
            }

            /* arc for second node set */
            if (azi[1][0] < azi[1][1]) {
                let az = azi[1][1] - D2R;
                while (az > azi[1][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp2[i] = x0 + radius_size * si;
                    i += 1;
                    yp2[i] = y0 + radius_size * co;
                    az -= D2R;
                }
            } else {
                let az = azi[1][1] + D2R;
                while (az < azi[1][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp2[i] = x0 + radius_size * si;
                    i += 1;
                    yp2[i] = y0 + radius_size * co;
                    az += D2R;
                }
            }

            collect.push(xy2patch(xp2.slice(0, i), yp2.slice(0, i), res, xy));
            colors.push(rgb1);
            return [colors, collect];
        }

        /* ---- n == 2 case ------------------------------------------------ */
        if (n === 2) {
            /* first concatenate x3 → xp1 */
            let i = 0;
            for (i = 0; i < j3; i++) {
                xp1[i] = x3[i];
                yp1[i] = y3[i];
            }
            i -= 1;

            /* then x → xp1 */
            for (let ii = 0; ii < j; ii++) {
                xp1[i] = x[ii];
                i++;
                yp1[i] = y[ii];
            }

            /* big_iso shortcut --------------------------------------------- */
            if (big_iso) {
                let ii = j2 - 1;
                while (ii >= 0) {
                    xp1[i] = x2[ii];
                    i += 1;
                    yp1[i] = y2[ii];
                    ii -= 1;
                }
                collect.push(xy2patch(xp1.slice(0, i), yp1.slice(0, i), res, xy));
                colors.push(rgb1);
                return [colors, collect];
            }

            /* unwrap azi[2][0] wrt azi[0][1] */
            if (azi[2][0] - azi[0][1] > Math.PI) {
                azi[2][0] -= Math.PI * 2.0;
            } else if (azi[0][1] - azi[2][0] > Math.PI) {
                azi[2][0] += Math.PI * 2.0;
            }

            /* arc between azi[0][1] and azi[2][0] */
            if (azi[2][0] < azi[0][1]) {
                let az = azi[0][1] - D2R;
                while (az > azi[2][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp1[i] = x0 + radius_size * si;
                    i += 1;
                    yp1[i] = y0 + radius_size * co;
                    az -= D2R;
                }
            } else {
                let az = azi[0][1] + D2R;
                while (az < azi[2][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp1[i] = x0 + radius_size * si;
                    i += 1;
                    yp1[i] = y0 + radius_size * co;
                    az += D2R;
                }
            }

            collect.push(xy2patch(xp1.slice(0, i), yp1.slice(0, i), res, xy));
            colors.push(rgb1);

            /* copy x2 → xp2 */
            for (i = 0; i < j2; i++) {
                xp2[i] = x2[i];
                yp2[i] = y2[i];
            }
            i -= 1;

            /* unwrap azi[1][0] wrt azi[1][1] */
            if (azi[1][0] - azi[1][1] > Math.PI) {
                azi[1][0] -= Math.PI * 2.0;
            } else if (azi[1][1] - azi[1][0] > Math.PI) {
                azi[1][0] += Math.PI * 2.0;
            }

            /* arc for second patch */
            if (azi[1][0] < azi[1][1]) {
                let az = azi[1][1] - D2R;
                while (az > azi[1][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp2[i] = x0 + radius_size * si;
                    i += 1;
                    yp2[i] = y0 + radius_size * co;
                    az -= D2R;
                }

            } else {
                let az = azi[1][1] + D2R;
                while (az < azi[1][0]) {
                    const si = Math.sin(az);
                    const co = Math.cos(az);
                    xp2[i] = x0 + radius_size * si;
                    i += 1;
                    yp2[i] = y0 + radius_size * co;
                    az += D2R;
                }
            }

            collect.push(xy2patch(xp2.slice(0, i), yp2.slice(0, i), res, xy));
            colors.push(rgb1);
            return [colors, collect];
        }

        /* ---- fallback (should never happen) ----------------------------- */
        return [colors, collect];

    }

    // import necessary classes and helper

    class BeachballComponent extends HTMLElement {
        static get observedAttributes() {
            return [
                'focal-mechanism',
                'width',
                'height',
                'tension-color',
                'bg-color',
                'edge-color',
                'line-width',
                'size'
            ];
        }

        constructor() {
            super();
            const shadow = this.attachShadow({ mode: 'open' });

            // styles for host & canvas
            const style = document.createElement('style');
            style.textContent = `
      :host {
        display: inline-block;
        position: relative;
        width: var(--bb-width, 200px);
        height: var(--bb-height, 200px);
      }
      canvas {
        width: 100%;
        height: 100%;
        display: block;
      }
    `;
            shadow.appendChild(style);

            // only the canvas
            this.canvas = document.createElement('canvas');
            shadow.appendChild(this.canvas);
        }

        attributeChangedCallback(name, oldVal, newVal) {
            switch (name) {
                case 'focal-mechanism':
                    try {
                        this.focalMechanism = JSON.parse(newVal);
                    } catch (err) {
                        console.error('⚠️ setter threw:', err);
                        console.log('raw attribute:', newVal);

                        // console.log(newVal);
                        // console.log(JSON.parse(newVal));
                        console.warn('Invalid JSON for focal-mechanism');
                    }
                    break;
                case 'width': {
                    const w = parseInt(newVal) || 200;
                    this.canvas.width = w;
                    if (!this.hasAttribute('height')) {
                        this.canvas.height = w;
                        this.style.setProperty('--bb-height', `${w}px`);
                    }
                    this.style.setProperty('--bb-width', `${w}px`);
                    break;
                }
                case 'height': {
                    const h = parseInt(newVal);
                    this.canvas.height = h;
                    this.style.setProperty('--bb-height', `${h}px`);
                    break;
                }
                case 'tension-color':
                    this.tensionColor = newVal;
                    break;
                case 'bg-color':
                    this.backgroundColor = newVal;
                    break;
                case 'edge-color':
                    this.edgeColor = newVal;
                    break;
                case 'line-width':
                    this.lineWidth = parseFloat(newVal);
                    break;
                case 'size':
                    this.size = Math.max(100, parseInt(newVal));
                    break;
            }
            this._draw();
        }

        set focalMechanism(fm) {
            let mt = null, np1 = null;
            if (Array.isArray(fm) && fm.length === 6) {
                mt = new MomentTensor(...fm, 0);
                np1 = mt2plane(mt);
            } else if (Array.isArray(fm) && fm.length === 3) {
                np1 = new NodalPlane(...fm);
            } else {
                throw new TypeError("Wrong input value for 'fm'.");
            }
            this._mt = mt;
            this._np1 = np1;
            this._draw();
        }
        get focalMechanism() {
            return this._mt ? this._mt : this._np1;
        }

        connectedCallback() {
            // defaults
            const w = this.hasAttribute('width') ? +this.getAttribute('width') : 200;
            this.canvas.width = w;
            this.style.setProperty('--bb-width', `${w}px`);

            const h = this.hasAttribute('height') ? +this.getAttribute('height') : w;
            this.canvas.height = h;
            this.style.setProperty('--bb-height', `${h}px`);

            this.tensionColor = this.getAttribute('tension-color') || '#87CEEB';
            this.backgroundColor = this.getAttribute('bg-color') || '#FFFFFF';
            this.edgeColor = this.getAttribute('edge-color') || '#000000';
            this.lineWidth = +this.getAttribute('line-width') || 2;
            this.size = Math.max(100, +this.getAttribute('size') || 100);

            this._draw();
        }

        _drawPatches(ctx, colors, patches, fill = true, drawLine = false) {
            ctx.save();

            // --- flip y-axis once for ALL patches ---
            // Adjust coordinate system: Matplotlib uses a bottom-left origin with y-axis up,
            // whereas Canvas uses a top-left origin with y-axis down. This translates the
            // origin to the bottom and flips the y-axis so patches render with correct symmetry.
            ctx.translate(0, ctx.canvas.height);
            ctx.scale(1, -1);

            ctx.lineWidth = this.lineWidth;
            ctx.strokeStyle = this.edgeColor;

            patches.forEach((patch, i) => {
                ctx.fillStyle = colors[i];
                ctx.beginPath();
                if (patch.type === 'ellipse') {
                    ctx.ellipse(
                        patch.xy[0], patch.xy[1],
                        patch.width / 2, patch.height / 2,
                        0, 0, 2 * Math.PI
                    );

                } else {
                    patch.vertices.forEach(([x, y], idx) =>
                        idx ? ctx.lineTo(x, y) : ctx.moveTo(x, y)
                    );
                }
                ctx.closePath();
                if (fill)
                    ctx.fill();
                if (drawLine)
                    ctx.stroke();
            });
            ctx.restore();
        }


        /**
         * Project an axis onto a beachball diagram.
         *
         * @param {number} strike – strike angle in degrees
         * @param {number} dip    – dip angle in degrees
         * @param {number} [R=90] – reference radius (pixels or arbitrary units)
         * @param {number} [margin=20] – extra margin added to the projected radius
         * @returns {[number, number]} [x, y] coordinates in the plotting plane
         */
        _axisOnBeachball(strike, dip, R = 90, margin = 20) {
            const az = strike * Math.PI / 180;              // radians
            const r = (90 - dip) * R / 90 + margin;        // projected radius
            const y = Math.cos(az) * r;
            const x = Math.sin(az) * r;
            return [x, y];                                  // same order as Python
        }


        /**
        * Return an array of fill colors based on the raw `colors` array:
        * if colors[i] === 'b' use tensionColor, else backgroundColor.
        */
        _getFillColors(colors, patches) {
            return patches.map((_, i) =>
                colors[i] === 'b' ? this.tensionColor : this.backgroundColor
            );
        }

        _draw() {
            if (!this._np1) return;
            const ctx = this.canvas.getContext('2d');

            // fill background
            ctx.fillStyle = this.backgroundColor;
            ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

            // prepare border style
            ctx.lineWidth = this.lineWidth;   // e.g. 2
            ctx.strokeStyle = this.edgeColor;    // e.g. 'black'

            // compute center and radius (or radii)
            const cx = this.canvas.width / 2;
            const cy = this.canvas.height / 2;

            // // draw the main border (circle or ellipse)
            // ctx.beginPath();
            // ctx.ellipse(cx, cy, cx - ctx.lineWidth, cy - ctx.lineWidth, 0, 0, 2 * Math.PI);
            // ctx.stroke();

            // now draw internal patches
            const xy = [cx, cy];
            const diameter = [this.canvas.width - 2 * ctx.lineWidth, this.canvas.height - 2 * ctx.lineWidth]; //r * 2;
            let colors, patches;

            if (this._mt) {
                const [t, n, p] = mt2axes(this._mt.normalized);

                // radius of the plotting circle in canvas units:
                const R = Math.min(this.canvas.width, this.canvas.height) / 2 - this.lineWidth;
                // same margin you used in Python:
                const MARGIN = 0; //18

                // 1) project to 2-D coordinates (matplotlib frame)
                const [xT, yT] = this._axisOnBeachball(t.strike, t.dip, R, MARGIN);
                const [xP, yP] = this._axisOnBeachball(p.strike, p.dip, R, MARGIN);

                if (Math.abs(n.val) < EPSILON && Math.abs(t.val + p.val) < EPSILON) {
                    // pure DC
                    [colors, patches] = plotDC(this._np1, this.size, xy, diameter);
                    const fillColors = this._getFillColors(colors, patches);
                    this._drawPatches(ctx, fillColors, patches, true, true);

                } else {
                    // MT + DC overlay
                    [colors, patches] = plotMT(t, n, p, this.size, true, 0, 0, xy, diameter);
                    const fillMT = this._getFillColors(colors, patches);
                    this._drawPatches(ctx, fillMT, patches, true, false);

                    [colors, patches] = plotDC(this._np1, this.size, xy, diameter);
                    const fillDC = this._getFillColors(colors, patches);
                    this._drawPatches(ctx, fillDC, patches, false, true);
                }

                // 2) flip to canvas coordinates (because we use y-down on canvas)
                const canvasXT = cx + xT;
                const canvasYT = cy - yT;   // minus sign flips y
                const canvasXP = cx + xP;
                const canvasYP = cy - yP;

                // 3) draw bold centred text
                ctx.save();
                ctx.font = 'bold 18px sans-serif';
                ctx.fillStyle = 'black';
                ctx.textAlign = 'center';
                ctx.textBaseline = 'middle';
                // “z-order” on canvas == call order, so draw last
                ctx.fillText('T', canvasXT, canvasYT);
                ctx.fillText('P', canvasXP, canvasYP);

                ctx.restore();

            } else {
                // fallback to DC only
                [colors, patches] = plotDC(this._np1, this.size, xy, diameter);
                const fillColors = this._getFillColors(colors, patches);
                this._drawPatches(ctx, fillColors, patches, true, true);
            }
        }


    }

    customElements.define('beachball-component', BeachballComponent);

}));
