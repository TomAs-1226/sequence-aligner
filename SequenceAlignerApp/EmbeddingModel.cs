using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Windows.AI.MachineLearning;

namespace SequenceAlignerApp;

/// <summary>
/// A transformer protein language model that turns a protein into a 192-number
/// "embedding". Runs with Windows ML on the GPU or NPU (CPU fallback). Trained
/// from scratch on ~16,000 real proteins. App-only feature.
/// </summary>
public sealed class EmbeddingModel
{
    private const string AA = "ACDEFGHIKLMNPQRSTVWY";
    private const int FIXED = 200;
    private readonly LearningModelSession? _session;

    public string Provider { get; private set; } = "none";
    public bool Ready => _session != null;

    public EmbeddingModel()
    {
        string path = Path.Combine(AppContext.BaseDirectory, "plm_embed.onnx");
        if (!File.Exists(path)) return;
        var model = LearningModel.LoadFromFilePath(path);
        foreach (var (kind, name) in new[]
                 {
                     (LearningModelDeviceKind.DirectXHighPerformance, "GPU / NPU"),
                     (LearningModelDeviceKind.Cpu, "CPU"),
                 })
        {
            try { _session = new LearningModelSession(model, new LearningModelDevice(kind)); Provider = name; break; }
            catch { }
        }
    }

    public float[]? Embed(string seq)
    {
        if (_session == null || seq.Length == 0) return null;
        var ids = new long[FIXED];
        var mask = new float[FIXED];
        for (int i = 0; i < FIXED; i++)
        {
            if (i < seq.Length)
            {
                int k = AA.IndexOf(char.ToUpperInvariant(seq[i]));
                ids[i] = k < 0 ? 22 : k;   // 22 = unknown
                mask[i] = 0f;
            }
            else { ids[i] = 20; mask[i] = 1f; }   // 20 = pad
        }
        var b = new LearningModelBinding(_session);
        b.Bind("seq", TensorInt64Bit.CreateFromArray(new long[] { 1, FIXED }, ids));
        b.Bind("mask", TensorFloat.CreateFromArray(new long[] { 1, FIXED }, mask));
        var r = _session.Evaluate(b, "0");
        var v = ((TensorFloat)r.Outputs["embedding"]).GetAsVectorView();
        var arr = new float[v.Count];
        for (int i = 0; i < v.Count; i++) arr[i] = v[i];
        return arr;
    }

    public static float Cosine(float[] a, float[] b)
    {
        float dot = 0, na = 0, nb = 0;
        for (int i = 0; i < a.Length; i++) { dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i]; }
        return dot / (MathF.Sqrt(na) * MathF.Sqrt(nb) + 1e-6f);
    }
}

/// <summary>A small built-in database of proteins with precomputed embeddings.</summary>
public sealed class ProteinDatabase
{
    public sealed record Entry(string Family, string Name, float[] Emb);
    public readonly List<Entry> Entries = new();

    public ProteinDatabase()
    {
        string path = Path.Combine(AppContext.BaseDirectory, "plm_db.tsv");
        if (!File.Exists(path)) return;
        foreach (var line in File.ReadLines(path))
        {
            var parts = line.Split('\t');
            if (parts.Length < 3) continue;
            var raw = parts[2].Split(',');
            var emb = new float[raw.Length];
            for (int i = 0; i < raw.Length; i++)
                emb[i] = float.Parse(raw[i], CultureInfo.InvariantCulture);
            Entries.Add(new Entry(parts[0], parts[1], emb));
        }
    }
}
