import { xy2patch } from '../src/beachball.js';

describe('xy2patch', () => {
    test('scalar res with translation', () => {
        const x = [0, 1, 1];
        const y = [0, 0, 1];
        const p = xy2patch(x, y, 1, [0, 0]);

        expect(p.vertices).toEqual([
            [0, 0],
            [1, 0],
            [1, 1],
        ]);
        expect(p.codes).toEqual(['M', 'L', 'Z']);
    });

    test('elliptical res scaling', () => {
        const x = [0, 1];
        const y = [0, 1];
        const p = xy2patch(x, y, [2, 3], [10, -5]);

        expect(p.vertices[1][0]).toBeCloseTo(1 * 2 + 10); // 12
        expect(p.vertices[1][1]).toBeCloseTo(1 * 3 - 5);  // -2
    });
});
