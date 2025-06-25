import { tdl } from '../src/beachball.js';

describe('tdl', () => {
    test('vector an=[1,0,0], bn=[0,1,0] → (270, 90, 180)', () => {
        const an = [1, 0, 0];
        const bn = [0, 1, 0];
        const [ft, fd, fl] = tdl(an, bn);
        expect(ft).toBeCloseTo(270);
        expect(fd).toBeCloseTo(90);
        expect(fl).toBeCloseTo(180);
    });
});
