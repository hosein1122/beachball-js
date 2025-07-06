
// mt2axes.js:

import { PrincipalAxis } from "../classes.js";
import { eigenSym3 } from "../math/eigenSym3.js";
import { R2D } from "../math/trigHelpers.js";

/**
 * Calculates principal axes T, N, P from a moment tensor.
 *
 * @param {MomentTensor} mt
 * @returns {[PrincipalAxis, PrincipalAxis, PrincipalAxis]}
 */
export function mt2axes(mt) {
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
