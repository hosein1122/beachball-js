// rollup.config.js
import resolve from '@rollup/plugin-node-resolve';
import { terser } from 'rollup-plugin-terser';
import commonjs from '@rollup/plugin-commonjs';

const terserOptions = {
  compress: {
    drop_console: true,
    passes: 2
  },
  mangle: {
    toplevel: true,
  },
  format: {
    comments: false,
  }
};

export default {
  // input: 'src/index.js',
  input: 'src/beachball-component.js',
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
      // plugins: [terser(terserOptions)]
      plugins: [terser(terserOptions)]
    },
    {
      file: 'dist/beachball.esm.js',
      format: 'es'        // ES Module bundle
    },
    {
      file: 'dist/beachball.esm.min.js',
      format: 'es',
      // plugins: [terser(terserOptions)]
      plugins: [terser(terserOptions)]
    }
  ],
  plugins: [
    resolve(),
    commonjs(),
  ]
};
