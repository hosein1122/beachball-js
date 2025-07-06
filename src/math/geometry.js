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
export function strikeDip(n, e, u) {
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
