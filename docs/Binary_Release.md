# Binary prerelease installation

Download the archive for your platform and `SHA256SUMS` from the [v0.1.0-rc.5 prerelease](https://github.com/zvanjak/MML-Packages/releases/tag/v0.1.0-rc.5). Verify the downloaded archive against its entry in `SHA256SUMS` before extracting it. On Linux, `sha256sum -c SHA256SUMS` checks all three archives when downloaded together; on macOS use `shasum -a 256 -c SHA256SUMS`. On Windows, compare `Get-FileHash` with the corresponding entry. Extract the archive and use the directory named `MML-Packages-0.1.0-rc.5-<platform>-<abi>` as the installation prefix.

| Archive suffix | Tested environment | ABI requirement |
| --- | --- | --- |
| `windows-x64-msvc143-md.zip` | Windows x64, MSVC 19.44 / v143 | Release, x64, dynamic MSVC runtime (`/MD`) |
| `linux-x64-glibc-libstdcxx.tar.gz` | Ubuntu 24.04 x86-64, GCC 13.3 | glibc 2.39; libstdc++ with `_GLIBCXX_USE_CXX11_ABI=1` |
| `macos-arm64-libcxx.tar.gz` | macOS 14 arm64, AppleClang 15 | macOS deployment target 14.0; libc++ |

These are tested build environments, not a promise of compatibility with older operating systems or different C++ standard libraries. All platforms use C++20, matching package and MML Core headers, and five static libraries. On Windows, do not mix the archive with `/MT` or a different architecture. The same headers and libraries must come from one extracted archive.

## CMake consumption

Point `CMAKE_PREFIX_PATH` at the extracted prefix and use `find_package(MMLPackages REQUIRED CONFIG)`. Link against `MML::Packages` for all five packages, or select `MML::Fourier`, `MML::Optimization`, `MML::PDE`, `MML::Statistics`, and `MML::Symbolic` individually. `MML::Core` and `MML::Ext` provide the bundled MML Core and extension headers. Installed targets propagate the static-build definition `MML_STATIC`; do not compile package implementation files into the consumer. The `docs_demos/` application also supports an extracted prefix through `MML_INSTALLED_PREFIX`.

The prefix contains `include/`, `lib/`, `lib/cmake/MMLPackages/`, `licenses/`, and `share/mml-packages/release-manifest.json`. The manifest records the source commit, vendored MML snapshot digest, compiler, architecture, runtime and archive version; check it when troubleshooting a mismatch. The package and bundled MML Core are MIT licensed; both required notices are in the archive's `licenses/` directory. Keep those notices with redistributed copies.

See the [release notes](Release_Notes_v0.1.0-rc.5.md) for exact provenance, checksums and known limitations.