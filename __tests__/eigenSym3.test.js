// __tests__/eigenSym3.test.js
import { eigenSym3 } from "../src/math/eigenSym3.js";

describe("eigenSym3", () => {
  test("zero matrix → zeros + identity", () => {
    const A = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    const { values, vectors } = eigenSym3(A);
    expect(values).toEqual([0, 0, 0]);
    expect(vectors).toEqual([
      [0, 0, 1],
      [0, 1, 0],
      [1, 0, 0],
    ]);
  });

  test("diagonal matrix → diag entries", () => {
    const A = [
      [5, 0, 0],
      [0, 2, 0],
      [0, 0, -1],
    ];
    const { values: desc, vectors: vecDesc } = eigenSym3(A, "desc");

    // desc order: max→mid→min
    expect(desc).toEqual([5, 2, -1]);
    // basis vectors (in any sign) — we'll check abs
    expect(vecDesc[0].map(Math.abs)).toEqual([1, 0, 0]);
    expect(vecDesc[1].map(Math.abs)).toEqual([0, 1, 0]);
    expect(vecDesc[2].map(Math.abs)).toEqual([0, 0, 1]);

    const { values: asc, vectors: vecAsc } = eigenSym3(A, "asc");

    expect(asc).toEqual([-1, 2, 5]);
  });

  test("random symmetric vs NumPy reference", () => {
    // small reproducible symmetric
    const A = [
      [4, 1, 2],
      [1, 5, 3],
      [2, 3, 6],
    ];
    const { values: valsJS, vectors: vecJS } = eigenSym3(A, "asc");

    let ref = [2.19439717, 3.38677016, 9.41883268];
    valsJS.forEach((v, i) => {
      expect(v).toBeCloseTo(ref[i], 6);
    });

    ref = [-0.4412247, -0.57735027, 0.68701342];
    vecJS[0].forEach((v, i) => {
      expect(v).toBeCloseTo(ref[i], 6);
    });

    ref = [0.81558342, -0.57735027, 0.03860509];
    vecJS[1].forEach((v, i) => {
      expect(v).toBeCloseTo(ref[i], 6);
    });

    ref = [0.37435872, 0.57735027, 0.7256185];
    vecJS[2].forEach((v, i) => {
      expect(v).toBeCloseTo(ref[i], 6);
    });

    // // reference from NumPy (precomputed)
    // const valsNP = [8.2469796, 4.0000000, 2.7530204]; // truncated
    // // compare within tolerance
    // valsJS.forEach((v, i) => {
    //   expect(v).toBeCloseTo(valsNP[i], 6);
    // });
    // // Check orthonormality: V·Vᵀ ≈ I
    // const V = vecJS;
    // for (let i = 0; i < 3; i++) {
    //   for (let j = 0; j < 3; j++) {
    //     const dot = V[i][0] * V[j][0]
    //       + V[i][1] * V[j][1]
    //       + V[i][2] * V[j][2];
    //     expect(dot).toBeCloseTo(i === j ? 1 : 0, 6);
    //   }
    // }
  });
});
