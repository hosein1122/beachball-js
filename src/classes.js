import { Matrix } from 'ml-matrix';

/* ──────────────────────────────────────────────────────────
 * PrincipalAxis
 * --------------------------------------------------------*/
export class PrincipalAxis {
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
export class NodalPlane {
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
export class MomentTensor {
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

    get mtMatrix() {
        return new Matrix(this.mt);
    }

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
