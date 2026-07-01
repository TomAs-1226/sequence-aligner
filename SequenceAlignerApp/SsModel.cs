using System;
using System.IO;
using System.Text;
using Windows.AI.MachineLearning;

namespace SequenceAlignerApp;

/// <summary>
/// A small neural network that predicts per-residue secondary structure
/// (helix / sheet / coil) from a protein sequence. It runs with Windows ML,
/// which uses the GPU or NPU when available and falls back to the processor.
/// The model was trained on real protein structures from the Protein Data Bank.
/// This is an app-only feature.
/// </summary>
public sealed class SsModel
{
    private const string AA = "ACDEFGHIKLMNPQRSTVWY";
    private readonly LearningModelSession? _session;

    public string Provider { get; private set; } = "none";
    public bool Ready => _session != null;

    public SsModel()
    {
        string path = Path.Combine(AppContext.BaseDirectory, "ss_model.onnx");
        if (!File.Exists(path)) return;
        var model = LearningModel.LoadFromFilePath(path);
        // prefer the GPU/NPU, fall back to the processor
        foreach (var (kind, name) in new[]
                 {
                     (LearningModelDeviceKind.DirectXHighPerformance, "GPU / NPU"),
                     (LearningModelDeviceKind.Cpu, "CPU"),
                 })
        {
            try
            {
                _session = new LearningModelSession(model, new LearningModelDevice(kind));
                Provider = name;
                break;
            }
            catch { }
        }
    }

    /// <summary>Per-residue prediction as a string of 'H' (helix), 'E' (sheet), 'C' (coil).</summary>
    public string Predict(string seq)
    {
        if (_session == null || seq.Length == 0) return "";
        var ids = new long[seq.Length];
        for (int i = 0; i < seq.Length; i++)
        {
            int idx = AA.IndexOf(char.ToUpperInvariant(seq[i]));
            ids[i] = idx < 0 ? AA.Length : idx;
        }
        var input = TensorInt64Bit.CreateFromArray(new long[] { 1, seq.Length }, ids);
        var binding = new LearningModelBinding(_session);
        binding.Bind("seq", input);
        var result = _session.Evaluate(binding, "0");
        var logits = (TensorFloat)result.Outputs["logits"];
        var v = logits.GetAsVectorView();   // flat [1 * L * 3], row-major

        var sb = new StringBuilder(seq.Length);
        for (int i = 0; i < seq.Length; i++)
        {
            int best = 0;
            float bv = v[i * 3];
            for (int c = 1; c < 3; c++)
                if (v[i * 3 + c] > bv) { bv = v[i * 3 + c]; best = c; }
            sb.Append("HEC"[best]);
        }
        return sb.ToString();
    }
}
