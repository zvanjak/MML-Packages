# MML-Packages v0.1.0-rc.5

First three-platform binary prerelease. Provides matching MML Core and package headers, five static package libraries, relocatable CMake targets, and both MIT license notices. Built and validated in the private release workflow on Windows x64, Ubuntu 24.04 x86-64, and macOS 14 arm64. All platform jobs and the aggregate checksum/provenance verifier passed, including full repository tests, extracted-prefix consumers and documentation demos.

**Source revision:** `9a446287269e4263fa1aa3842b745026c08e4796` (clean checkout on all three runners). **Vendored MML snapshot SHA-256:** `b0e5272595d190bca938e6c62a9f137321cdce4cc3ff366ecb476e205303781b` (content digest recorded in each archive's release manifest). The vendored MML snapshot is not yet tied to a formal upstream 2.0 tag.

| Asset | SHA-256 |
| --- | --- |
| `MML-Packages-0.1.0-rc.5-linux-x64-glibc-libstdcxx.tar.gz` | `f4f0d4064eefa69bc353983101328cb327df6f78d8d9ffa4fabb6b38fea50267` |
| `MML-Packages-0.1.0-rc.5-macos-arm64-libcxx.tar.gz` | `5fb6244225ee1dc69be400732387c1b1dabcf00f2a71ede310210deac9fbb40c` |
| `MML-Packages-0.1.0-rc.5-windows-x64-msvc143-md.zip` | `57be592b2d812d1337b941a405c5f3f8621afdac12c7fe8433580fde604e1f71` |

Download `SHA256SUMS` alongside the three archives to verify their bytes. See the [installation and compatibility guide](Binary_Release.md) before linking these static libraries.

**Known limitation:** Independently rebuilt Windows `.lib` files are not yet proven byte-for-byte identical, even though re-packaging the same staged prefix is deterministic. The archives above passed their integrity, ABI, and extracted-consumer checks; do not interpret these checksums as a cross-build reproducibility guarantee.