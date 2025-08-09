# SimpleFlac

A single-file, dependency-free FLAC decoder written in modern C# (.NET 8+). The core of this repository is `FlacDecoder.cs` which can be easily integrated into other projects.

A console demo is included, but the main purpose of this repository is to share the `FlacDecoder` implementation.

## Features
- Pure C# implementation with no `unsafe` code or native dependencies.
- Decoder is entirely self-contained in a single file (`FlacDecoder.cs`).
- Optimized for performance while maintaining readability.
- Permissively licensed (MIT).

## Limitations
- 8 bits per sample is not supported.
- Bit depths that are not a multiple of 8 (e.g. 20 bits per sample) are not supported.

## Usage

### Library
```csharp
using FlacDecoder decoder = new("input.flac");
while (decoder.DecodeFrame()) {
    // Access decoded data via:
    // * BufferSamples: long[0..ChannelCount][0..BufferSampleCount]
    // * BufferBytes: byte[0..BufferByteCount]
}
```

### Console Demo
```sh
dotnet run -c Release -- input.flac output.wav
```

## Credits
- C# port and enhancements: J.D. Purcell
- Original Java implementation: Project Nayuki ([nayuki.io](https://www.nayuki.io/page/simple-flac-implementation))
