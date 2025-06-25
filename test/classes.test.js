import { PrincipalAxis, NodalPlane, MomentTensor } from '../src/classes.js';

describe('PrincipalAxis', () => {
    test('basic props', () => {
        const a = new PrincipalAxis(1.3, 20, 50);
        expect(a.val).toBeCloseTo(1.3);
        expect(a.strike).toBeCloseTo(20);
        expect(a.dip).toBeCloseTo(50);
    });
});

describe('NodalPlane', () => {
    test('basic props', () => {
        const p = new NodalPlane(13, 20, 50);
        expect(p.strike).toBe(13);
        expect(p.dip).toBe(20);
        expect(p.rake).toBe(50);
    });
});

describe('MomentTensor', () => {
    test('6-component constructor', () => {
        const mt = new MomentTensor([1, 1, 0, 0, 0, -1], 26);
        expect(mt.xx).toBe(1);
        expect(mt.yz).toBe(-1);
        expect(mt.expo).toBe(26);
    });

    test('full-matrix constructor', () => {
        const m = [
            [1, 0, 0],
            [0, 1, -1],
            [0, -1, 0],
        ];
        const mt = new MomentTensor(m, 26);
        expect(mt.yy).toBe(1);
        expect(mt.xz).toBe(0);
    });

    test('normalized tensor', () => {
        const mt = new MomentTensor([1, 1, 1, 0, 0, 0], 0).normalized;
        const norm = Math.sqrt(mt.xx ** 2 + mt.yy ** 2 + mt.zz ** 2);
        expect(norm).toBeCloseTo(1);
    });
});
