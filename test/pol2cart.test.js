// test/pol2cart.test.js
import { pol2cart } from '../src/beachball.js';

describe('pol2cart', () => {
    test('scalar inputs', () => {
        const [x, y] = pol2cart(Math.PI / 4, 1);
        expect(x).toBeCloseTo(Math.SQRT1_2);
        expect(y).toBeCloseTo(Math.SQRT1_2);
    });

    test('array inputs', () => {
        const th = [0, Math.PI / 2, Math.PI];
        const r = [1, 1, 1];
        const [x, y] = pol2cart(th, r);

        // x should be [1, 0, -1] within precision
        expect(x[0]).toBeCloseTo(1);
        expect(x[1]).toBeCloseTo(0);
        expect(x[2]).toBeCloseTo(-1);

        // y should be [0, 1, 0] within precision
        expect(y[0]).toBeCloseTo(0);
        expect(y[1]).toBeCloseTo(1);
        expect(y[2]).toBeCloseTo(0);
    });
});
