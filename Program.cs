using System;
using System.Diagnostics;
using System.IO;

namespace SimpleFlac;

static class Program {
	static int Main(string[] args) {
		if (args.Length != 2) {
			Console.WriteLine("Usage: SimpleFlac input.flac output.wav");
			return 1;
		}

		string inputPath = args[0];
		string outputPath = args[1];
		Console.WriteLine($"Decoding: {Path.GetFileName(inputPath)}");
		long startTimestamp = Stopwatch.GetTimestamp();
		TimeSpan audioDuration = DecodeFile(inputPath, outputPath);
		TimeSpan elapsedTime = Stopwatch.GetElapsedTime(startTimestamp);
		Console.WriteLine($"Finished in {elapsedTime.TotalSeconds:F3} seconds ({audioDuration / elapsedTime:F1}x realtime)");
		return 0;
	}

	static TimeSpan DecodeFile(string inputPath, string outputPath) {
		using FlacDecoder decoder = new(inputPath);
		long streamSampleCount = decoder.StreamSampleCount ??
			throw new NotSupportedException("Stream sample count cannot be unknown.");
		using WavWriter writer = new(outputPath, decoder.SampleRate, decoder.ChannelCount, decoder.BitsPerSample, streamSampleCount);
		while (decoder.DecodeFrame()) {
			writer.WriteData(decoder.BufferBytes.AsSpan(0, decoder.BufferByteCount));
		}
		writer.EndData();
		return TimeSpan.FromSeconds((double)streamSampleCount / decoder.SampleRate);
	}
}
