// mt2plane.js :

import { NodalPlane } from "../classes.js";
import { eigenSym3 } from "../math/eigenSym3.js";
import { tdl } from "../math/geometry.js";

/**
 * Computes a NodalPlane object from a given MomentTensor.
 *
 * @param {MomentTensor} mt - Moment tensor instance
 * @returns {NodalPlane}
 */
export function mt2plane(mt) {
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
