// beachball.js
// Core functions converted from beachball.py
// ───────────────────────────────────────────
import { MomentTensor, NodalPlane, PrincipalAxis } from './classes.js';
import { Matrix, EigenvalueDecomposition } from 'ml-matrix';
import { D2R, R2D, EPSILON } from './constants.js';

// ───────────────────────────────────────────
// 1. pol2cart  ← ported from beachball.py::pol2cart
// --------------------------------------------------------
/**
* Notes:
* - If either input is an array, the other is broadcast to the same length.
* - Throws if array lengths mismatch.
 *
 * @param {number | number[]} th angle(s) in radians
 * @param {number | number[]} r radius/radii
 * @returns {[number | number[], number | number[]]} (same shape as inputs)
 */
export function pol2cart(th, r) {
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

// ───────────────────────────────────────────
// 2. strikeDip  ← ported from beachball.py::strike_dip
// --------------------------------------------------------
/**
 * Finds strike and dip of a plane given its normal-vector components (n, e, u).
 * Mirrors the exact logic of the original Python implementation.
 *
 * @param {number} n - North component of the normal vector
 * @param {number} e - East  component of the normal vector
 * @param {number} u - Up    component of the normal vector
 * @returns {[number, number]} [strikeDeg, dipDeg] in degrees
 */
export function strikeDip(n, e, u) {
    const r2d = 180 / Math.PI;

    // Ensure upward-pointing normal
    if (u < 0) {
        n = -n;
        e = -e;
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

// ───────────────────────────────────────────
// 3. xy2patch  ← port of beachball.py::xy2patch
// --------------------------------------------------------
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
export function xy2patch(x, y, res, xy) {
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








    // // Normalize res into [sx, sy]
    // let sx, sy;
    // if (Array.isArray(res)) {
    //     if (res.length !== 2) throw new Error('res must be scalar or length-2');
    //     [sx, sy] = res;
    // } else {
    //     sx = sy = res;
    // }

    // // Transform points
    // const vertices = x.map((xi, idx) => {
    //     const yi = y[idx];
    //     return [xi * sx + xy[0], yi * sy + xy[1]];
    // });

    // // Build drawing codes: 'M' + ('L' * (n-2)) + 'Z'
    // const codes = ['M', ...Array(Math.max(vertices.length - 2, 0)).fill('L'), 'Z'];

    // return { vertices, codes };
}

// ───────────────────────────────────────────
// 4. auxPlane  → port of beachball.py::aux_plane
// --------------------------------------------------------
/**
 * Calculates strike, dip, and rake of the auxiliary (second) plane given the
 * first plane's strike (s1), dip (d1), and rake (r1). All angles in degrees.
 *
 * @param {number} s1  Strike of first plane (deg)
 * @param {number} d1  Dip    of first plane (deg)
 * @param {number} r1  Rake   of first plane (deg)
 * @returns {[number, number, number]} [strike, dip, rake] of second plane (deg)
 */
export function auxPlane(s1, d1, r1) {
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

// ───────────────────────────────────────────
// 5. tdl  → port of beachball.py::tdl
//--------------------------------------------------------
/**
 * Calculates three angles (ft, fd, fl) from vectors an and bn.
 * Used internally by mt2plane. All angles in degrees.
 *
 * @param {number[]} an - vector a (length 3)
 * @param {number[]} bn - vector b (length 3)
 * @returns {[number, number, number]} [ft, fd, fl] in degrees
 */
export function tdl(an, bn) {
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

// ───────────────────────────────────────────
// 6. mt2plane  → port of beachball.py::mt2plane
//--------------------------------------------------------
/**
 * Computes a NodalPlane object from a given MomentTensor.
 *
 * @param {MomentTensor} mt - Moment tensor instance
 * @returns {NodalPlane}
 */
export function mt2plane(mt) {
    // Eigen decomposition
    const eig = new EigenvalueDecomposition(new Matrix(mt.mt));
    let values = eig.realEigenvalues;
    let vectors = eig.eigenvectorMatrix.to2DArray();

    // Re-order to match Python indexing trick [d1,d0,d2]
    values = [values[1], values[0], values[2]];
    vectors = [
        [vectors[1][1], -vectors[1][0], -vectors[1][2]],
        [vectors[2][1], -vectors[2][0], -vectors[2][2]],
        [-vectors[0][1], vectors[0][0], vectors[0][2]],
    ];

    // Identify max/min eigenvalue indices
    const imax = values.indexOf(Math.max(...values));
    const imin = values.indexOf(Math.min(...values));

    // Build ae & an vectors
    const ae = [
        (vectors[0][imax] + vectors[0][imin]) / Math.sqrt(2),
        (vectors[1][imax] + vectors[1][imin]) / Math.sqrt(2),
        (vectors[2][imax] + vectors[2][imin]) / Math.sqrt(2),
    ];
    const an = [
        (vectors[0][imax] - vectors[0][imin]) / Math.sqrt(2),
        (vectors[1][imax] - vectors[1][imin]) / Math.sqrt(2),
        (vectors[2][imax] - vectors[2][imin]) / Math.sqrt(2),
    ];

    // Normalize
    const norm = v => Math.hypot(...v);
    const aer = norm(ae);
    const anr = norm(an);

    const aeN = ae.map(x => x / aer);
    let anN;
    if (anr === 0) {
        anN = [Number.NaN, Number.NaN, Number.NaN];
    } else {
        anN = an.map(x => x / anr);
    }

    // Ensure downward component negative
    const an1 = anN[2] <= 0 ? anN : anN.map(x => -x);
    const ae1 = anN[2] <= 0 ? aeN : aeN.map(x => -x);

    // Angles via tdl
    const [ft, fd, fl] = tdl(an1, ae1);

    return new NodalPlane(360 - ft, fd, 180 - fl);
}

// ───────────────────────────────────────────
// 7. mt2axes  → port of beachball.py::mt2axes
//--------------------------------------------------------
/**
 * Calculates principal axes T, N, P from a moment tensor.
 *
 * @param {MomentTensor} mt
 * @returns {[PrincipalAxis, PrincipalAxis, PrincipalAxis]}
 */
export function mt2axes(mt) {
    // Symmetric eigen-decomposition
    const evd = new EigenvalueDecomposition(new Matrix(mt.mt));
    const d = evd.realEigenvalues;          // eigenvalues
    const v = evd.eigenvectorMatrix.to2DArray(); // eigenvectors (columns)

    // Convert arrays for convenience: v[row][col]
    // pl, az calculations
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
        if (pl[i] <= 0) {
            pl[i] = -pl[i];
            az[i] += Math.PI;
        }
        if (az[i] < 0) az[i] += 2 * Math.PI;
        if (az[i] > 2 * Math.PI) az[i] -= 2 * Math.PI;
    }

    pl = pl.map(rad => rad * R2D);
    az = az.map(rad => rad * R2D);

    const t = new PrincipalAxis(d[2], az[2], pl[2]);
    const n = new PrincipalAxis(d[1], az[1], pl[1]);
    const p = new PrincipalAxis(d[0], az[0], pl[0]);

    return [t, n, p];
}


/**
 * Convert Strike-Dip-Rake (degrees)  →  6-component moment tensor
 * in the Up-South-East (r,θ,φ) CMT convention.
 *
 * @returns [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
 */
export function sdr2mt(strikeDeg, dipDeg, rakeDeg, M0 = 1) {
    const toRad = Math.PI / 180;

    const phi = strikeDeg * toRad;   // ϕ
    const delta = dipDeg * toRad;  // δ
    const lambda = rakeDeg * toRad;  // λ

    // handy trig
    const sinD = Math.sin(delta);
    const cosD = Math.cos(delta);
    const sin2D = Math.sin(2 * delta);
    const sinL = Math.sin(lambda);
    const cosL = Math.cos(lambda);
    const sin2L = Math.sin(2 * lambda);
    const cos2L = Math.cos(2 * lambda);
    const sinP = Math.sin(phi);
    const cosP = Math.cos(phi);
    const sin2P = Math.sin(2 * phi);
    const cos2P = Math.cos(2 * phi);

    // components (Aki-Richards formulation, CMT convention)
    const Mrr = -sinD * sin2L;
    const Mtt = sin2D * sinL * sin2P + sinD * sinD * cos2L * (1 - cos2P);
    const Mpp = -sin2D * sinL * sin2P - sinD * sinD * cos2L * (1 + cos2P);
    const Mrt = sinD * cosL * cosP + cosD * sinL * sinP;
    const Mrp = sinD * cosL * sinP - cosD * sinL * cosP;
    const Mtp = -sin2D * cosL + sinD * sinD * sin2L * cos2P;

    return [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp].map(v => v * M0);
}

// Example:
// console.log(sdr2mt(240, 46, -119));



// ───────────────────────────────────────────
// 8. plotDC → JavaScript port of plot_dc (double-couple beachball)
// ---------------------------------------------------------------
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
export function plotDC(np1, size = 200, xy = [0, 0], width = 200) {
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
    const [s_2, d_2, _r_2] = auxPlane(s_1, d_1, r_1);

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

export function plotMT(T, N, P, size = 200, plot_zerotrace = true, x0 = 0, y0 = 0, xy = [0, 0], width = 200) {

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

// ───────────────────────────────────────────
// 10. beach → high-level port of beachball.py::beach
//--------------------------------------------------------
/**
 * Build a generic, renderer-agnostic “collection” for a beach-ball plot.
 * No direct drawing here—just data and style metadata.
 *
 * @param {Array|MomentTensor|NodalPlane} fm
 *   Focal mechanism (3- or 6-element array, or MomentTensor/NodalPlane instance).
 * @param {Object} [opts]
 * @param {number} [opts.linewidth=2]
 * @param {string} [opts.facecolor='b']
 * @param {string} [opts.bgcolor='w']
 * @param {string} [opts.edgecolor='k']
 * @param {number} [opts.alpha=1.0]
 * @param {Array<number>} [opts.xy=[0,0]]
 * @param {number|Array<number>} [opts.width=200]
 * @param {number} [opts.size=100]
 * @param {boolean} [opts.nofill=false]
 * @param {number} [opts.zorder=100]
 * @param {*} [opts.axes=null]   // kept for API compatibility, unused here
 *
 * @returns {Object}
 *   • patches: Array of data-only patch descriptors (ellipse or polygon)  
 *   • fillColors: Array of fill colors (or null if nofill)  
 *   • edgecolor, linewidth, alpha, zorder  
 *   • xy, width, size  
 *   • plotDcUsed: boolean flag which algorithm was used
 */
export function beach(fm, {
    // linewidth = 2,//--
    facecolor = 'b',
    bgcolor = 'w',
    // edgecolor = 'k',//--
    // alpha = 1.0,//--
    xy = [0, 0],
    width = 200,
    size = 100,
    nofill = false,
    // zorder = 100//--
} = {}) {
    // ─── 1) Ensure minimum resolution ───
    if (size < 100) size = 100;

    // ─── 2) Normalize width to [wx, wy] ───
    if (!Array.isArray(width) || width.length !== 2) {
        width = [width, width];
    }

    // ─── 3) Parse fm into np1 and optional MomentTensor mt ───
    let mt = null, np1 = null;
    if (fm instanceof MomentTensor) {
        mt = fm;
        np1 = mt2plane(mt);
    } else if (fm instanceof NodalPlane) {
        np1 = fm;
    } else if (Array.isArray(fm) && fm.length === 6) {
        mt = new MomentTensor(...fm, 0);
        np1 = mt2plane(mt);
    } else if (Array.isArray(fm) && fm.length === 3) {
        np1 = new NodalPlane(...fm);
    } else {
        throw new TypeError("Wrong input value for 'fm'.");
    }

    // ─── 4) Choose DC or MT algorithm ───
    let colors, patches, plotDcUsed = true;
    if (mt) {
        const [t, n, p] = mt2axes(mt.normalized);
        if (Math.abs(n.val) < EPSILON && Math.abs(t.val + p.val) < EPSILON) {
            [colors, patches] = plotDC(np1, size, xy, width);
        } else {
            [colors, patches] = plotMT(t, n, p, size,
                true, 0, 0, xy, width);

            plotDcUsed = false;
        }
    } else {
        [colors, patches] = plotDC(np1, size, xy, width);
    }



    // ─── 5) Build fillColors array ───
    const fillColors = patches.map((_, i) => {
        if (nofill)
            return null;
        return colors[i] === 'b' ? facecolor : bgcolor;
    });

    // ─── 6) Return data-only collection ───
    return {
        patches,        // [{ type:'ellipse', xy:[x,y], width, height } | { vertices:[[x,y],…], res:[rx,ry], center:[cx,cy] }]
        fillColors,     // [color|string|null]
        // edgecolor,      // string
        // linewidth,      // number
        // alpha,          // number
        // zorder,         // number
        xy,             // [x,y]
        width,          // [wx,wy]
        size,           // number
        plotDcUsed      // boolean
    };
}


