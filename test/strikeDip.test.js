import { strikeDip } from '../src/beachball.js';

describe('strikeDip', () => {
    test('normal vector [1, 0, 1] → strike 270°, dip 45°', () => {
        const [strike, dip] = strikeDip(1, 0, 1);
        expect(strike).toBeCloseTo(270);
        expect(dip).toBeCloseTo(45);
    });

    test('flipped vector (u < 0) handled correctly', () => {
        const [strike, dip] = strikeDip(-1, 0, -1); // same plane, inverted
        expect(strike).toBeCloseTo(270);
        expect(dip).toBeCloseTo(45);
    });
});
