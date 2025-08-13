using System;
using System.IO;

#nullable enable

namespace SimpleFlac;

public class WavWriter : IDisposable {
	private readonly BinaryWriter _output;
	private readonly int _sampleRate;
	private readonly int _channelCount;
	private readonly int _bitsPerSample;
	private readonly int _blockAlign;
	private readonly long _sampleCount;
	private bool _writePaddingByte;
	private long _writtenSampleCount;
	private bool _ended;

	public WavWriter(string path, int sampleRate, int channelCount, int bitsPerSample, long sampleCount) {
		_sampleRate = sampleRate;
		_channelCount = channelCount;
		_bitsPerSample = bitsPerSample;
		_sampleCount = sampleCount;
		_blockAlign = channelCount * ((bitsPerSample + 7) / 8);
		_output = new BinaryWriter(new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 65536, FileOptions.SequentialScan));
		try {
			WriteHeaders();
		}
		catch {
			_output.Dispose();
			throw;
		}
	}

	private void WriteHeaders() {
		long dataSize = _sampleCount * _blockAlign;
		_writePaddingByte = dataSize % 2 != 0;
		long fileSize = 44 + dataSize + (_writePaddingByte ? 1 : 0);
		if (fileSize > UInt32.MaxValue)
			throw new NotSupportedException("File size exceeds the maximum.");
		_output.Write("RIFF"u8);
		_output.Write((uint)(fileSize - 8));
		_output.Write("WAVE"u8);
		_output.Write("fmt "u8);
		_output.Write((uint)16); // fmt chunk size
		_output.Write((ushort)1); // PCM
		_output.Write((ushort)_channelCount);
		_output.Write((uint)_sampleRate);
		_output.Write((uint)(_sampleRate * _blockAlign));
		_output.Write((ushort)_blockAlign);
		_output.Write((ushort)_bitsPerSample);
		_output.Write("data"u8);
		_output.Write((uint)dataSize);
	}

	public void WriteData(ReadOnlySpan<byte> data) {
		if (_ended)
			throw new InvalidOperationException("Cannot write data after ending.");
		if (data.Length % _blockAlign != 0)
			throw new ArgumentException("Data length must be a multiple of block align.");
		int dataSampleCount = data.Length / _blockAlign;
		if (_writtenSampleCount + dataSampleCount > _sampleCount)
			throw new InvalidDataException("Too many samples written.");
		_output.Write(data);
		_writtenSampleCount += dataSampleCount;
	}

	public void EndData() {
		if (_ended)
			throw new InvalidOperationException("Data was already ended.");
		if (_writtenSampleCount < _sampleCount)
			throw new InvalidDataException("Not enough samples written.");
		if (_writePaddingByte) {
			_output.Write((byte)0);
		}
		_ended = true;
	}

	public void Dispose() {
		_output.Dispose();
	}
}
