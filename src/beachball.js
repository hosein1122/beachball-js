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
    // Normalize res into [sx, sy]
    let sx, sy;
    if (Array.isArray(res)) {
        if (res.length !== 2) throw new Error('res must be scalar or length-2');
        [sx, sy] = res;
    } else {
        sx = sy = res;
    }

    // Transform points
    const vertices = x.map((xi, idx) => {
        const yi = y[idx];
        return [xi * sx + xy[0], yi * sy + xy[1]];
    });

    // Build drawing codes: 'M' + ('L' * (n-2)) + 'Z'
    const codes = ['M', ...Array(Math.max(vertices.length - 2, 0)).fill('L'), 'Z'];

    return { vertices, codes };
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

// ───────────────────────────────────────────
// 8. plotDC → JavaScript port of plot_dc (double-couple beachball)
// ---------------------------------------------------------------
/**
 * Builds pressure & tension polygons for a double-couple beachball.
 *
 * @param {NodalPlane} plane            – first nodal plane
 * @param {object}     opts
 * @param {number}     opts.radius      – X-radius in px (mandatory)
 * @param {number=}    opts.radiusY     – Y-radius; defaults to radius
 * @param {number[]}   opts.center      – [cx, cy] center of ball
 * @param {number=}    opts.phiStep     – angular step in radians (default 0.01)
 * @returns {{pressure:Patch, tension:Patch}}
 *
 * Where Patch = { vertices:number[][], codes:string[] }
 */
export function plotDC(plane, {
    radius,
    radiusY = radius,
    center = [0, 0],
    phiStep = 0.01
} = {}) {
    const s1 = plane.strike;
    const d1 = plane.dip;
    let r1 = plane.rake;

    let m = 0;
    if (r1 > 180) { r1 -= 180; m = 1; }
    if (r1 < 0) { r1 += 180; m = 1; }

    const [s2, d2] = auxPlane(s1, d1, r1);

    const d1c = Math.min(d1, 89.9999);
    const d2c = Math.min(d2, 89.9999);

    // build phi array
    const phi = [];
    for (let p = 0; p < Math.PI; p += phiStep) phi.push(p);

    // length functions
    const lCalc = dipDeg => {
        const A = 90 - dipDeg;
        return phi.map(p =>
            Math.sqrt(
                (A ** 2) /
                (Math.sin(p) ** 2 + (Math.cos(p) ** 2) * (A ** 2) / (90 ** 2))
            )
        );
    };

    const l1 = lCalc(d1c);
    const l2 = lCalc(d2c);

    const patches = [];

    for (const m_ of [((m + 1) % 2), m]) {
        // first arc
        const [x1Arr, y1Arr] = pol2cart(
            phi.map(p => p + s1 * D2R),
            l1
        );

        // build th1 & th2
        let th1 = [], th2 = [], x2Arr, y2Arr, inc = 1;

        if (m_ === 1) {
            const lo = s1 - 180, hi = s2;
            if (lo > hi) inc = -1;
            for (let t = lo; inc > 0 ? t < hi : t > hi; t += inc) th1.push(t);
            [x2Arr, y2Arr] = pol2cart(
                phi.map(p => p + s2 * D2R),
                l2
            );
            for (let t = s2 + 180; inc > 0 ? t > s1 : t < s1; t -= inc) th2.push(t);
        } else {
            const hi = s1 - 180, lo = s2 - 180;
            if (lo > hi) inc = -1;
            for (let t = hi; inc > 0 ? t < lo : t > lo; t -= inc) th1.push(t);
            [x2Arr, y2Arr] = pol2cart(
                phi.map(p => p + s2 * D2R),
                l2
            );
            x2Arr.reverse(); y2Arr.reverse();
            for (let t = s2; inc > 0 ? t < s1 : t > s1; t += inc) th2.push(t);
        }

        // correct radii arrays to match angles
        const ones1 = new Array(th1.length).fill(90);
        const ones2 = new Array(th2.length).fill(90);

        const [xs1, ys1] = pol2cart(
            th1.map(t => t * D2R),
            ones1
        );
        const [xs2, ys2] = pol2cart(
            th2.map(t => t * D2R),
            ones2
        );

        // concatenate
        const x = [...x1Arr, ...xs1, ...x2Arr, ...xs2];
        const y = [...y1Arr, ...ys1, ...y2Arr, ...ys2];

        // scale to radius
        const vx = x.map(v => v * radius / 90);
        const vy = y.map(v => v * radiusY / 90);

        // build patch
        patches.push(xy2patch(vy, vx, [1, 1], center));
    }

    return { tension: patches[0], pressure: patches[1] };
}

// ───────────────────────────────────────────
// 9. plotMT → port of beachball.py::plot_mt
//--------------------------------------------------------
/**
 * Builds patches for a moment-tensor–based beachball.
 *
 * @param {PrincipalAxis} T
 * @param {PrincipalAxis} N
 * @param {PrincipalAxis} P
 * @param {object} opts
 * @param {number} opts.radius      – radius in px (controls overall size)
 * @param {number=} opts.radiusY    – y-axis radius; defaults to radius
 * @param {number[]} opts.center    – [cx, cy] center offset
 * @param {boolean=} opts.zeroTrace – if true, force zero isotropic part (default true)
 * @returns {{tension:{vertices:number[][],codes:string[]}, pressure:{vertices:number[][],codes:string[]}}}
 */
export function plotMT(T, N, P, {
    radius,
    radiusY = radius,
    center = [0, 0],
    zeroTrace = true
} = {}) {
    // unpack
    const a = [T.strike, N.strike, P.strike];
    const p = [T.dip, N.dip, P.dip];
    const v = [T.val, N.val, P.val];

    // isotropic part
    let vi = (v[0] + v[1] + v[2]) / 3;
    for (let i = 0; i < 3; i++) v[i] -= vi;

    // early return iso-only
    const mag = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if (Math.abs(mag) < EPSILON) {
        return {
            tension: { vertices: [], codes: [] },
            pressure: { vertices: [], codes: [] }
        };
    }

    // choose principal axes indices
    const d = Math.abs(v[0]) >= Math.abs(v[2]) ? 0 : 2;
    const m = d === 0 ? 2 : 0;
    if (zeroTrace) vi = 0;
    const f = -v[1] / v[d];
    const iso = vi / v[d];

    // early return no nodes
    if (iso < -1 || iso > 1 - f) {
        return {
            tension: { vertices: [], codes: [] },
            pressure: { vertices: [], codes: [] }
        };
    }

    // precompute trig
    const spd = Math.sin(p[d] * D2R),
        cpd = Math.cos(p[d] * D2R),
        spb = Math.sin(p[1] * D2R),
        cpb = Math.cos(p[1] * D2R),
        spm = Math.sin(p[m] * D2R),
        cpm = Math.cos(p[m] * D2R),
        sad = Math.sin(a[d] * D2R),
        cad = Math.cos(a[d] * D2R),
        sab = Math.sin(a[1] * D2R),
        cab = Math.cos(a[1] * D2R),
        sam = Math.sin(a[m] * D2R),
        cam = Math.cos(a[m] * D2R);

    // allocate arrays
    const azi = Array.from({ length: 3 }, () => [0, 0]);
    const x = new Array(400).fill(0),
        y = new Array(400).fill(0),
        x2 = new Array(400).fill(0),
        y2 = new Array(400).fill(0),
        x3 = new Array(400).fill(0),
        y3 = new Array(400).fill(0),
        xp1 = new Array(800).fill(0),
        yp1 = new Array(800).fill(0),
        xp2 = new Array(400).fill(0),
        yp2 = new Array(400).fill(0);

    let bigIso = 0, j = 1, j2 = 0, j3 = 0, nIdx = 0, azp = 0;

    // loop 0..359
    for (let i = 0; i < 360; i++) {
        const fir = i * D2R;
        const s2a = (2 + 2 * iso) / (3 + (1 - 2 * f) * Math.cos(2 * fir));
        if (s2a > 1) { bigIso++; continue; }
        const al = Math.asin(Math.sqrt(s2a)),
            sfi = Math.sin(fir),
            cfi = Math.cos(fir),
            san = Math.sin(al),
            can = Math.cos(al);
        const xz = can * spd + san * sfi * spb + san * cfi * spm;
        let xn = can * cpd * cad + san * sfi * cpb * cab + san * cfi * cpm * cam,
            xe = can * cpd * sad + san * sfi * cpb * sab + san * cfi * cpm * sam;
        let takeoff, az;
        if (Math.abs(xn) < EPSILON && Math.abs(xe) < EPSILON) {
            takeoff = 0; az = 0;
        } else {
            az = Math.atan2(xe, xn);
            if (az < 0) az += 2 * Math.PI;
            takeoff = Math.acos(xz / Math.hypot(xz, xn, xe));
        }
        if (takeoff > Math.PI / 2) {
            takeoff = Math.PI - takeoff;
            az += Math.PI;
            if (az > 2 * Math.PI) az -= 2 * Math.PI;
        }
        const r = Math.SQRT2 * Math.sin(takeoff / 2),
            si = Math.sin(az),
            co = Math.cos(az);

        if (i === 0) {
            azi[0][0] = az;
            x[0] = center[0] + radius * r * si;
            y[0] = center[1] + radius * r * co;
            azp = az;
        } else {
            const dAz = Math.abs(Math.abs(az - azp) - Math.PI);
            if (dAz < 10 * D2R) {
                azi[nIdx][1] = azp;
                nIdx++;
                azi[nIdx][0] = az;
            }
            if (nIdx === 0) {
                x[j] = center[0] + radius * r * si;
                y[j] = center[1] + radius * r * co;
                j++;
            } else if (nIdx === 1) {
                xp1[j2] = center[0] + radius * r * si;
                yp1[j2] = center[1] + radius * r * co;
                j2++;
            } else {
                xp2[j3] = center[0] + radius * r * si;
                yp2[j3] = center[1] + radius * r * co;
                j3++;
            }
            azp = az;
        }
    }
    azi[nIdx][1] = azp;

    // prepare patches
    const patches = [];

    // central circle
    patches.push(
        xy2patch(
            [center[1], center[1]],
            [center[0], center[0]],
            [radiusY, radius],
            [0, 0]
        )
    );

    // main fill
    if (nIdx === 0) {
        patches.push(xy2patch(x.slice(0, j), y.slice(0, j), [1, 1], [0, 0]));
    } else if (nIdx === 1) {
        // combine xp1 (len j2) then segment from x/j, then xp2
        const seg = xp1.slice(0, j2).concat(x.slice(0, j)).concat(x2.slice(0, j3));
        const segY = yp1.slice(0, j2).concat(y.slice(0, j)).concat(y2.slice(0, j3));
        patches.push(xy2patch(seg, segY, [1, 1], [0, 0]));
    } else {
        // nIdx===2
        // combine xp2, xp1, x
        const seg = xp2.slice(0, j3)
            .concat(xp1.slice(0, j2))
            .concat(x.slice(0, j));
        const segY = yp2.slice(0, j3)
            .concat(yp1.slice(0, j2))
            .concat(y.slice(0, j));
        patches.push(xy2patch(seg, segY, [1, 1], [0, 0]));
    }

    // assign colors: tension first, pressure next
    const rgb1 = v[1] < 0 ? 'b' : 'w';
    const rgb2 = v[1] < 0 ? 'w' : 'b';
    const [tPatch, pPatch] = patches;
    const makeCodes = verts => {
        if (verts.length === 0) return [];
        const c = ['M'];
        for (let i = 1; i < verts.length - 1; i++) c.push('L');
        c.push('Z');
        return c;
    };

    return {
        tension: { vertices: tPatch.vertices, codes: makeCodes(tPatch.vertices) },
        pressure: { vertices: pPatch.vertices, codes: makeCodes(pPatch.vertices) }
    };
}

// ───────────────────────────────────────────
// 10. beach → high-level port of beachball.py::beach
//--------------------------------------------------------
/**
 * Produces two styled patches (tension & pressure) for a beachball.
 *
 * @param {MomentTensor|NodalPlane|number[]} fm
 *   - MomentTensor instance OR
 *   - NodalPlane instance OR
 *   - Array[3] ([strike,dip,rake]) OR
 *   - Array[6] (moment tensor components)
 * @param {object} opts
 * @param {string} opts.facecolor  – color for tension quadrants (default 'b')
 * @param {string} opts.bgcolor    – color for pressure quadrants (default 'w')
 * @param {string} opts.edgecolor  – stroke color for outlines (default 'k')
 * @param {number} opts.alpha      – fill opacity [0..1] (default 1.0)
 * @param {number[]} opts.center   – [cx, cy] in px (default [0,0])
 * @param {number} opts.radius     – radius in px (mandatory)
 * @param {number=} opts.radiusY   – y-axis radius; defaults to radius
 * @param {number=} opts.phiStep   – angular step for DC (default 0.01)
 * @param {boolean=} opts.zeroTrace– force zero isotropic for MT (default true)
 * @returns {{
 *   tension:  {vertices:number[][], codes:string[], facecolor:string, edgecolor:string, alpha:number},
 *   pressure: {vertices:number[][], codes:string[], facecolor:string, edgecolor:string, alpha:number}
 * }}
 */
export function beach(fm, {
    facecolor = 'b',
    bgcolor = 'w',
    edgecolor = 'k',
    alpha = 1.0,
    center = [0, 0],
    radius,
    radiusY,
    phiStep = 0.01,
    zeroTrace = true
} = {}) {
    if (radius == null) {
        throw new Error('beach: opts.radius is required');
    }
    if (radiusY == null) radiusY = radius;

    let plane = null, mt = null;

    // determine input type
    if (fm instanceof MomentTensor) {
        mt = fm;
        plane = mt2plane(mt);
    } else if (fm instanceof NodalPlane) {
        plane = fm;
    } else if (Array.isArray(fm) && fm.length === 6) {
        mt = new MomentTensor(fm[0], fm[1], fm[2], fm[3], fm[4], fm[5], 0);
        plane = mt2plane(mt);
    } else if (Array.isArray(fm) && fm.length === 3) {
        plane = new NodalPlane(fm[0], fm[1], fm[2]);
    } else {
        throw new TypeError("beach: 'fm' must be MomentTensor, NodalPlane, or array of length 3 or 6");
    }

    let tensionPatch, pressurePatch;

    if (mt) {
        // check if pure double-couple or full moment-tensor
        const [T, N, P] = mt2axes(mt.normalized);
        if (Math.abs(N.val) < EPSILON && Math.abs(T.val + P.val) < EPSILON) {
            // pure double-couple from nodal plane
            ({ tension: tensionPatch, pressure: pressurePatch } =
                plotDC(plane, { radius, radiusY, center, phiStep }));
        } else {
            // use full moment-tensor logic
            ({ tension: tensionPatch, pressure: pressurePatch } =
                plotMT(T, N, P, { radius, radiusY, center, zeroTrace }));
        }
    } else {
        // only nodal-plane → always DC
        ({ tension: tensionPatch, pressure: pressurePatch } =
            plotDC(plane, { radius, radiusY, center, phiStep }));
    }

    // helper to map 'b'/'w' to actual colors
    const mapCol = c => (c === 'b' ? facecolor : bgcolor);

    // style the two patches
    const stylePatch = (patch, cDummy) => ({
        vertices: patch.vertices,
        codes: patch.codes,
        facecolor: mapCol(cDummy),
        edgecolor,
        alpha
    });

    return {
        tension: stylePatch(tensionPatch, 'b'),
        pressure: stylePatch(pressurePatch, 'w')
    };
}
