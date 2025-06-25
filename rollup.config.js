// rollup.config.js
import resolve from '@rollup/plugin-node-resolve';
import { terser } from 'rollup-plugin-terser';
import commonjs from '@rollup/plugin-commonjs';

export default {
  input: 'src/index.js',
  output: [
    {
      file: 'dist/beachball.umd.js',
      format: 'umd',
      name: 'Beachball'
    },
    {
      file: 'dist/beachball.umd.min.js',
      format: 'umd',
      name: 'Beachball',
      plugins: [terser()]
    },
    {
      file: 'dist/beachball.esm.js',
      format: 'es'        // ES Module bundle
    },
    {
      file: 'dist/beachball.esm.min.js',
      format: 'es',
      plugins: [terser()]
    }
  ],
  plugins: [
    resolve(),
    commonjs(),
  ]
};
