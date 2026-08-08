# Third-Party Notices

Glycine includes a locally patched copy of edlib_rs 0.1.2 under
vendor/edlib_rs-0.1.2.

- edlib_rs is distributed under the MIT or Apache-2.0 license.
- The bundled Edlib C/C++ implementation is distributed under the MIT license.
- Its upstream license is retained at
  vendor/edlib_rs-0.1.2/edlib-c/LICENSE.

The local patch raises the CMake compatibility declaration and adds lib64 to
the native-library search paths. These changes allow the dependency to build
with newer CMake versions and on Linux systems that install Edlib under
lib64.

Upstream projects:

- edlib-rs: https://github.com/jean-pierreBoth/edlib-rs
- Edlib: https://github.com/Martinsos/edlib
