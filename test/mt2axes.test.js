import { MomentTensor } from '../src/classes.js';
import { mt2axes } from '../src/beachball.js';

describe('mt2axes', () => {
    test('principal axes ordering and angles finite', () => {
        const mt = new MomentTensor([1, 1, -2, 0, 0, 0], 0);
        const [t, n, p] = mt2axes(mt);

        // Eigenvalue order: t >= n >= p (no absolute value)
        expect(t.val).toBeGreaterThanOrEqual(n.val);
        expect(n.val).toBeGreaterThanOrEqual(p.val);

        // Angles in valid ranges
        [t, n, p].forEach(ax => {
            expect(ax.strike).toBeGreaterThanOrEqual(0);
            expect(ax.strike).toBeLessThanOrEqual(360);
            expect(ax.dip).toBeGreaterThanOrEqual(0);
            expect(ax.dip).toBeLessThanOrEqual(90);
        });
    });
});
