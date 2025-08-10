// License: MIT
//
// Copyright (c) J.D. Purcell (C# port and enhancements)
// Copyright (c) Project Nayuki (Simple FLAC decoder in Java)
// https://www.nayuki.io/page/simple-flac-implementation
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

using System;
using System.Buffers.Binary;
using System.IO;
using System.Numerics;
using System.Security.Cryptography;
using System.Threading.Tasks;

#nullable enable

namespace SimpleFlac;

public class FlacDecoder : IDisposable {
	private readonly Options _options;
	private readonly BitReader _reader;
	private readonly IncrementalHash? _outputHasher;
	private Task _outputHasherTask = Task.CompletedTask;
	private byte[] _expectedOutputHash = new byte[16];

	public long? StreamSampleCount { get; private set; }
	public int SampleRate { get; private set; }
	public int ChannelCount { get; private set; }
	public int BitsPerSample { get; private set; }
	public int MaxSamplesPerFrame { get; private set; }

	public long[][] BufferSamples { get; private set; } = [[]];
	public byte[] BufferBytes { get; private set; } = [];
	public int BufferSampleCount { get; private set; }
	public int BufferByteCount { get; private set; }
	public long RunningSampleCount { get; private set; }

	public int BytesPerSample => BitsPerSample / 8;
	public int BlockAlign => BytesPerSample * ChannelCount;

	public FlacDecoder(Stream input, Options? options = null) {
		_options = options ?? new Options();
		_reader = new BitReader(input);
		try {
			ValidateOptions();
			ReadMetadata();
		}
		catch {
			_reader.Dispose();
			throw;
		}
		_outputHasher = _options.ValidateOutputHash ? IncrementalHash.CreateHash(HashAlgorithmName.MD5) : null;
	}

	public FlacDecoder(string path, Options? options = null)
		: this(new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 65536, FileOptions.SequentialScan), options)
	{
	}

	public void Dispose() {
		_reader.Dispose();
		if (_options.ValidateOutputHash) {
			_outputHasherTask.Wait();
			_outputHasher!.Dispose();
		}
	}

	private void ValidateOptions() {
		if (_options.ValidateOutputHash && !_options.ConvertOutputToBytes)
			throw new ArgumentException("Output hash validation requires conversion to bytes.");
	}

	private void ReadMetadata() {
		if (_reader.Read(32) != 0x664C6143)
			throw new InvalidDataException("FLAC stream marker not found.");

		bool foundLastMetadataBlock;
		do {
			foundLastMetadataBlock = _reader.Read(1) != 0;
			int type = (int)_reader.Read(7);
			int length = (int)_reader.Read(24);
			if (type == 0) {
				ReadStreaminfoBlock();
			}
			else {
				// Skip other blocks
				for (int i = 0; i < length; i++) {
					_reader.Skip(8);
				}
			}
		}
		while (!foundLastMetadataBlock);

		if (BufferSamples is null)
			throw new InvalidDataException("Stream info metadata block not found.");
	}

	private void ReadStreaminfoBlock() {
		_reader.Skip(16); // Minimum block size (samples)
		MaxSamplesPerFrame = (int)_reader.Read(16);
		_reader.Skip(24); // Minimum frame size (bytes)
		_reader.Skip(24); // Maximum frame size (bytes)
		SampleRate = (int)_reader.Read(20);
		ChannelCount = (int)_reader.Read(3) + 1;
		BitsPerSample = (int)_reader.Read(5) + 1;
		long streamSampleCount = (long)_reader.Read(36);
		for (int i = 0; i < 16; i++) {
			_expectedOutputHash[i] = (byte)_reader.Read(8);
		}

		if (BitsPerSample % 8 != 0)
			throw new NotSupportedException("Unsupported bit depth.");

		StreamSampleCount = streamSampleCount != 0 ? streamSampleCount : null;
		BufferSamples = new long[ChannelCount][];
		for (int ch = 0; ch < ChannelCount; ch++) {
			BufferSamples[ch] = new long[MaxSamplesPerFrame];
		}
		if (_options.ConvertOutputToBytes) {
			BufferBytes = new byte[MaxSamplesPerFrame * BlockAlign];
		}
	}

	public bool DecodeFrame() {
		if (_reader.HasReachedEnd) {
			if (StreamSampleCount is not null && RunningSampleCount != StreamSampleCount) {
				throw new InvalidDataException("Stream sample count is incorrect.");
			}

			if (_options.ValidateOutputHash) {
				Span<byte> actualHash = stackalloc byte[16];
				_outputHasherTask.Wait();
				_outputHasher!.GetCurrentHash(actualHash);

				if (!actualHash.SequenceEqual(_expectedOutputHash))
					throw new InvalidDataException("Output hash is incorrect.");
			}

			return false;
		}

		if (_reader.Read(15) != 0x7FFC)
			throw new InvalidDataException("Invalid frame sync code.");

		_reader.Skip(1); // Variable block size flag
		int blockSizeCode = (int)_reader.Read(4);
		int sampleRateCode = (int)_reader.Read(4);
		int channelLayout = (int)_reader.Read(4);
		int bitDepthCode = (int)_reader.Read(3);
		_reader.Skip(1); // Reserved bit

		// Coded number (sample or frame number)
		int codedNumberLeadingOnes = BitOperations.LeadingZeroCount(~(_reader.Read(8) << 56));
		for (int i = 1; i < codedNumberLeadingOnes; i++) {
			_reader.Skip(8);
		}

		int frameSampleCount = blockSizeCode switch {
			1 => 192,
			>= 2 and <= 5 => 576 << (blockSizeCode - 2),
			6 => (int)_reader.Read(8) + 1,
			7 => (int)_reader.Read(16) + 1,
			>= 8 and <= 15 => 256 << (blockSizeCode - 8),
			_ => throw new InvalidDataException("Reserved block size.")
		};

		int frameSampleRate = sampleRateCode switch {
			0 => SampleRate,
			>= 1 and <= 11 => SampleRateCodes[sampleRateCode],
			12 => (int)_reader.Read(8) * 1000,
			13 => (int)_reader.Read(16),
			14 => (int)_reader.Read(16) * 10,
			_ => throw new InvalidDataException("Reserved sample rate.")
		};

		int frameBitsPerSample = bitDepthCode switch {
			0 => BitsPerSample,
			>= 1 and <= 2 => 8 + ((bitDepthCode - 1) * 4),
			>= 4 and <= 6 => 16 + ((bitDepthCode - 4) * 4),
			7 => 32,
			_ => throw new InvalidDataException("Reserved bit depth.")
		};

		int frameChannelCount = channelLayout switch {
			>= 0 and <= 7 => channelLayout + 1,
			>= 8 and <= 10 => 2,
			_ => throw new InvalidDataException("Reserved channel layout.")
		};

		if (frameSampleCount > MaxSamplesPerFrame)
			throw new InvalidDataException("Frame sample count exceeds maximum.");

		if (frameSampleRate != SampleRate || frameBitsPerSample != BitsPerSample || frameChannelCount != ChannelCount)
			throw new NotSupportedException("Unsupported audio property change.");

		_reader.Skip(8); // Frame header CRC

		BufferSampleCount = frameSampleCount;
		RunningSampleCount += frameSampleCount;
		DecodeSubframes(_reader, BitsPerSample, channelLayout, BufferSamples, BufferSampleCount);
		_reader.AlignToByte();
		_reader.Skip(16); // Whole frame CRC

		if (_options.ConvertOutputToBytes) {
			_outputHasherTask.Wait();
			ConvertOutputToBytes(BitsPerSample, ChannelCount, BufferSamples, BufferSampleCount, BufferBytes);
			BufferByteCount = BufferSampleCount * BlockAlign;
		}

		if (_options.ValidateOutputHash) {
			_outputHasherTask = Task.Run(() => {
				_outputHasher!.AppendData(BufferBytes.AsSpan(0, BufferByteCount));
			});
		}

		return true;
	}

	private static void DecodeSubframes(BitReader reader, int bitsPerSample, int channelLayout, long[][] result, int blockSize) {
		if (channelLayout >= 0 && channelLayout <= 7) {
			for (int ch = 0; ch < result.Length; ch++) {
				DecodeSubframe(reader, bitsPerSample, result[ch].AsSpan(0, blockSize));
			}
		}
		else if (channelLayout >= 8 && channelLayout <= 10) {
			DecodeSubframe(reader, bitsPerSample + (channelLayout == 9 ? 1 : 0), result[0].AsSpan(0, blockSize));
			DecodeSubframe(reader, bitsPerSample + (channelLayout == 9 ? 0 : 1), result[1].AsSpan(0, blockSize));
			if (channelLayout == 8) {
				for (int i = 0; i < blockSize; i++) {
					result[1][i] = result[0][i] - result[1][i];
				}
			}
			else if (channelLayout == 9) {
				for (int i = 0; i < blockSize; i++) {
					result[0][i] += result[1][i];
				}
			}
			else if (channelLayout == 10) {
				for (int i = 0; i < blockSize; i++) {
					long side = result[1][i];
					long right = result[0][i] - (side >> 1);
					result[1][i] = right;
					result[0][i] = right + side;
				}
			}
		}
		else {
			throw new ArgumentOutOfRangeException(nameof(channelLayout));
		}
	}

	private static void DecodeSubframe(BitReader reader, int bitsPerSample, Span<long> result) {
		if (reader.Read(1) != 0)
			throw new InvalidDataException("Invalid subframe padding.");

		int type = (int)reader.Read(6);
		int shift = (int)reader.Read(1);
		if (shift == 1) {
			while (reader.Read(1) == 0) {
				shift++;
			}
		}
		bitsPerSample -= shift;

		if (type == 0) { // Constant coding
			long v = reader.ReadSigned(bitsPerSample);
			for (int i = 0; i < result.Length; i++) {
				result[i] = v;
			}
		}
		else if (type == 1) { // Verbatim coding
			for (int i = 0; i < result.Length; i++) {
				result[i] = reader.ReadSigned(bitsPerSample);
			}
		}
		else if (type >= 8 && type <= 12) {
			DecodeFixedPredictionSubframe(reader, type - 8, bitsPerSample, result);
		}
		else if (type >= 32 && type <= 63) {
			DecodeLinearPredictiveCodingSubframe(reader, type - 31, bitsPerSample, result);
		}
		else {
			throw new InvalidDataException("Reserved subframe type.");
		}

		if (shift != 0) {
			for (int i = 0; i < result.Length; i++) {
				result[i] <<= shift;
			}
		}
	}

	private static void DecodeFixedPredictionSubframe(BitReader reader, int predOrder, int bitsPerSample, Span<long> result) {
		for (int i = 0; i < predOrder; i++) {
			result[i] = reader.ReadSigned(bitsPerSample);
		}
		DecodeResiduals(reader, predOrder, result);
		if (predOrder != 0) {
			RestoreLinearPrediction(result, FixedPredictionCoefficients[predOrder], 0);
		}
	}

	private static void DecodeLinearPredictiveCodingSubframe(BitReader reader, int lpcOrder, int bitsPerSample, Span<long> result) {
		for (int i = 0; i < lpcOrder; i++) {
			result[i] = reader.ReadSigned(bitsPerSample);
		}
		int precision = (int)reader.Read(4) + 1;
		int shift = (int)reader.ReadSigned(5);
		Span<long> coefs = stackalloc long[lpcOrder];
		for (int i = coefs.Length - 1; i >= 0; i--) {
			coefs[i] = reader.ReadSigned(precision);
		}
		DecodeResiduals(reader, lpcOrder, result);
		RestoreLinearPrediction(result, coefs, shift);
	}

	private static void DecodeResiduals(BitReader reader, int warmup, Span<long> result) {
		int method = (int)reader.Read(2);
		if (method >= 2)
			throw new InvalidDataException("Reserved residual coding method.");
		int paramBits = method == 0 ? 4 : 5;
		int escapeParam = method == 0 ? 15 : 31;

		int partitionOrder = (int)reader.Read(4);
		int numPartitions = 1 << partitionOrder;
		if (result.Length % numPartitions != 0)
			throw new InvalidDataException("Block size not divisible by number of Rice partitions.");
		int partitionSize = result.Length / numPartitions;

		for (int i = 0; i < numPartitions; i++) {
			int start = i * partitionSize + (i == 0 ? warmup : 0);
			int end = (i + 1) * partitionSize;

			int param = (int)reader.Read(paramBits);
			if (param != escapeParam) {
				for (int j = start; j < end; j++) {
					result[j] = DecodeRice(reader, param);
				}
			}
			else {
				int numBits = (int)reader.Read(5);
				for (int j = start; j < end; j++) {
					result[j] = numBits != 0 ? reader.ReadSigned(numBits) : 0;
				}
			}
		}
	}

	private static void RestoreLinearPrediction(Span<long> result, ReadOnlySpan<long> coefs, int shift) {
		for (int i = 0; i < result.Length - coefs.Length; i++) {
			long sum = 0;
			for (int j = 0; j < coefs.Length; j++) {
				sum += result[i + j] * coefs[j];
			}
			result[i + coefs.Length] += sum >> shift;
		}
	}

	private static long DecodeRice(BitReader reader, int k) {
		ulong data = reader.RawBuffer;
		int leadingZeroCount = BitOperations.LeadingZeroCount(data);
		int quotientBitCount = leadingZeroCount + 1;
		int fullBitCount = quotientBitCount + k;
		if (fullBitCount > BitReader.BitsAvailableWorstCase) {
			return DecodeRiceFallback(reader, k);
		}
		ulong v = (ulong)leadingZeroCount << k;
		if (k != 0) {
			v |= (data << quotientBitCount) >> (64 - k);
		}
		reader.Skip(fullBitCount);
		// Apply sign from LSB
		return (int)(v >> 1) ^ -(int)(v & 1);
	}

	private static long DecodeRiceFallback(BitReader reader, int k) {
		int leadingZeroCount = 0;
		while (reader.Read(1) == 0) {
			leadingZeroCount++;
		}
		ulong v = (ulong)leadingZeroCount << k;
		if (k != 0) {
			v |= reader.Read(k);
		}
		// Apply sign from LSB
		return (int)(v >> 1) ^ -(int)(v & 1);
	}

	private static void ConvertOutputToBytes(int bitsPerSample, int channelCount, long[][] samples, int sampleCount, byte[] bytes) {
		int bytesPerSample = bitsPerSample / 8;
		int blockAlign = bytesPerSample * channelCount;
		for (int ch = 0; ch < channelCount; ch++) {
			long[] src = samples[ch];
			int offset = ch * bytesPerSample;
			if (bitsPerSample == 16) {
				Span<byte> byteSpan = bytes.AsSpan();
				for (int i = 0; i < sampleCount; i++) {
					BinaryPrimitives.WriteInt16LittleEndian(byteSpan.Slice(offset, 2), (short)src[i]);
					offset += blockAlign;
				}
			}
			else if (bitsPerSample == 24) {
				for (int i = 0; i < sampleCount; i++) {
					long s = src[i];
					bytes[offset    ] = (byte)s;
					bytes[offset + 1] = (byte)(s >> 8);
					bytes[offset + 2] = (byte)(s >> 16);
					offset += blockAlign;
				}
			}
			else if (bitsPerSample == 32) {
				Span<byte> byteSpan = bytes.AsSpan();
				for (int i = 0; i < sampleCount; i++) {
					BinaryPrimitives.WriteInt32LittleEndian(byteSpan.Slice(offset, 4), (int)src[i]);
					offset += blockAlign;
				}
			}
			else {
				// 8-bit is intentionally not supported because MD5 hash validation expects
				// signed data whereas a consumer of this class would probably want unsigned
				// data (e.g. to output to a WAV file).
				throw new NotSupportedException("Unsupported bit depth.");
			}
		}
	}

	private static readonly int[] SampleRateCodes = [
		0, 88200, 176400, 192000, 8000, 16000, 22050, 24000, 32000, 44100, 48000, 96000
	];

	private static readonly long[][] FixedPredictionCoefficients = [
		[],
		[1],
		[-1, 2],
		[1, -3, 3],
		[-1, 4, -6, 4]
	];

	private class BitReader : IDisposable {
		// Buffer replenish logic ensures that a full byte isn't missing
		public const int BitsAvailableWorstCase = 57;

		private readonly Stream _stream;
		private ulong _buffer;
		private int _bufferDeficitBits;
		private int _streamOverreadBytes;

		public BitReader(Stream stream) {
			_stream = stream;
			_bufferDeficitBits = 64;
			ReplenishBuffer();
		}

		public void Dispose() {
			_stream.Dispose();
		}

		public bool HasReachedEnd =>
			_streamOverreadBytes >= 8;

		public ulong RawBuffer =>
			_buffer;

		private void ReplenishBuffer() {
			while (_bufferDeficitBits >= 8) {
				int b = _stream.ReadByte();
				if (b == -1) {
					_streamOverreadBytes++;
					if (HasReachedEnd) {
						if (_bufferDeficitBits == 8) {
							// End was exactly reached; leave deficit so subsequent reads will throw
							return;
						}
						throw new EndOfStreamException();
					}
				}
				else {
					_buffer |= (ulong)b << (_bufferDeficitBits - 8);
				}
				_bufferDeficitBits -= 8;
			}
		}

		public void Skip(int numBits) {
			if (numBits < 1 || numBits > BitsAvailableWorstCase)
				throw new ArgumentOutOfRangeException(nameof(numBits));
			_buffer <<= numBits;
			_bufferDeficitBits += numBits;
			ReplenishBuffer();
		}

		public ulong Read(int numBits) {
			ulong x = _buffer >> (64 - numBits);
			Skip(numBits);
			return x;
		}

		public long ReadSigned(int numBits) {
			ulong x = Read(numBits);
			int shift = 64 - numBits;
			return (long)(x << shift) >> shift;
		}

		public void AlignToByte() {
			if (_bufferDeficitBits != 0) {
				Skip(8 - _bufferDeficitBits);
			}
		}
	}

	public class Options {
		public bool ConvertOutputToBytes { get; set; } = true;
		public bool ValidateOutputHash { get; set; } = true;
	}
}
