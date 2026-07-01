using System;
using System.Collections.Generic;
using System.Text;
using System.Threading.Tasks;
using Microsoft.UI.Xaml;
using Microsoft.UI.Xaml.Controls;
using Microsoft.UI.Xaml.Documents;
using Microsoft.UI.Xaml.Media;
using Microsoft.UI.Xaml.Shapes;
using Windows.ApplicationModel.DataTransfer;
using Windows.Storage;
using Windows.Storage.Pickers;
using WinRT.Interop;

namespace SequenceAlignerApp;

public sealed partial class MainPage : Page
{
    private static readonly SolidColorBrush GreenB = new(Windows.UI.Color.FromArgb(255, 0x2E, 0x7D, 0x46));
    private static readonly SolidColorBrush RedB = new(Windows.UI.Color.FromArgb(255, 0xC4, 0x3D, 0x2E));
    private static readonly SolidColorBrush GrayB = new(Windows.UI.Color.FromArgb(255, 0x8A, 0x8A, 0x8A));
    private const int MaxViewColumns = 1500;

    private bool _ready;
    private int _generation;
    private string _plainText = "";
    private bool _protein;

    public MainPage()
    {
        InitializeComponent();
        TypeBox.SelectedIndex = 0;
        ModeBox.SelectedIndex = 0;
        MatrixBox.SelectedIndex = 0;
        GapModelBox.SelectedIndex = 0;
        SampleBox.SelectedIndex = 0;
        Seq1Box.Text = SampleData.ToyDna1;
        Seq2Box.Text = SampleData.ToyDna2;
        _ready = true;
        UpdateVisibility();
        RunAlign();
    }

    // ---- input cleaning ----
    private static string Clean(string s)
    {
        var sb = new StringBuilder();
        foreach (var line in s.Split('\n'))
        {
            var t = line.Trim();
            if (t.Length == 0 || t.StartsWith(">")) continue;
            foreach (var ch in t)
                if (!char.IsWhiteSpace(ch)) sb.Append(char.ToUpperInvariant(ch));
        }
        return sb.ToString();
    }

    private void UpdateVisibility()
    {
        AffinePanel.Visibility = GapModelBox.SelectedIndex == 1 ? Visibility.Visible : Visibility.Collapsed;
        MatrixBox.IsEnabled = TypeBox.SelectedIndex == 1;
        bool dna = TypeBox.SelectedIndex == 0;
        MatchSlider.IsEnabled = dna;
        MismatchSlider.IsEnabled = dna;
    }

    // ---- event handlers ----
    private void Params_Changed(object sender, SelectionChangedEventArgs e)
    {
        if (!_ready) return;
        UpdateVisibility();
        RunAlign();
    }

    private void Slider_Changed(object sender, Microsoft.UI.Xaml.Controls.Primitives.RangeBaseValueChangedEventArgs e)
    {
        if (!_ready) return;
        MatchLabel.Text = $"Match reward: {MatchSlider.Value:0.#}";
        MismatchLabel.Text = $"Mismatch penalty: {MismatchSlider.Value:0.#}";
        GapLabel.Text = $"Gap penalty: {GapSlider.Value:0.#}";
        OpenLabel.Text = $"Gap open: {OpenSlider.Value:0.#}";
        ExtendLabel.Text = $"Gap extend: {ExtendSlider.Value:0.#}";
        RunAlign();
    }

    private void Align_Click(object sender, RoutedEventArgs e) => RunAlign();

    private void Sample_Changed(object sender, SelectionChangedEventArgs e)
    {
        if (!_ready) return;
        switch (SampleBox.SelectedIndex)
        {
            case 1: Seq1Box.Text = SampleData.ToyDna1; Seq2Box.Text = SampleData.ToyDna2; TypeBox.SelectedIndex = 0; break;
            case 2: Seq1Box.Text = SampleData.HbbHuman; Seq2Box.Text = SampleData.HbbGorilla; TypeBox.SelectedIndex = 1; break;
            case 3: Seq1Box.Text = SampleData.HbbHuman; Seq2Box.Text = SampleData.HbbZebrafish; TypeBox.SelectedIndex = 1; break;
            case 4: Seq1Box.Text = SampleData.SarsCovSpikeDna; Seq2Box.Text = SampleData.SarsCov2SpikeDna; TypeBox.SelectedIndex = 0; break;
            case 5: Seq1Box.Text = Aligner.Translate(SampleData.SarsCovSpikeDna); Seq2Box.Text = Aligner.Translate(SampleData.SarsCov2SpikeDna); TypeBox.SelectedIndex = 1; break;
            case 6: Seq1Box.Text = SampleData.SarsCovGenome; Seq2Box.Text = SampleData.SarsCov2Genome; TypeBox.SelectedIndex = 0; break;
            default: return;
        }
        UpdateVisibility();
        RunAlign();
    }

    private void Translate_Click(object sender, RoutedEventArgs e)
    {
        Seq1Box.Text = Aligner.Translate(Clean(Seq1Box.Text));
        Seq2Box.Text = Aligner.Translate(Clean(Seq2Box.Text));
        TypeBox.SelectedIndex = 1;
        UpdateVisibility();
        RunAlign();
    }

    private void Swap_Click(object sender, RoutedEventArgs e)
    {
        (Seq1Box.Text, Seq2Box.Text) = (Seq2Box.Text, Seq1Box.Text);
        RunAlign();
    }

    private void RevComp_Click(object sender, RoutedEventArgs e)
    {
        Seq2Box.Text = Aligner.ReverseComplement(Clean(Seq2Box.Text));
        RunAlign();
    }

    private async void Load3D_Click(object sender, RoutedEventArgs e)
    {
        ThreeDStatus.Text = "Working out which residues changed...";
        // find the changed residues from the two spike proteins, using our own aligner
        string p1 = Aligner.Translate(SampleData.SarsCovSpikeDna);
        string p2 = Aligner.Translate(SampleData.SarsCov2SpikeDna);
        var aln = await Task.Run(() => Aligner.Global(p1, p2, Aligner.FromMatrix(Aligner.Blosum62), 10));
        var changed = new List<int>();
        int resi = 0;
        for (int i = 0; i < aln.Row1.Length; i++)
        {
            char a = aln.Row1[i], b = aln.Row2[i];
            if (b == '-') continue;
            resi++;
            if (a == '-' || a != b) changed.Add(resi);
        }
        try
        {
            await Viewer3D.EnsureCoreWebView2Async();
            Viewer3D.NavigateToString(Build3DHtml(string.Join(",", changed)));
            ThreeDStatus.Text = $"{changed.Count} changed residues shown in red. Drag to rotate.";
        }
        catch
        {
            ThreeDStatus.Text = "Could not start the 3D viewer. The WebView2 runtime may be missing on this PC.";
        }
    }

    private static string Build3DHtml(string resiList)
    {
        return "<!DOCTYPE html><html><head><meta charset='utf-8'>" +
               "<script src='https://3Dmol.org/build/3Dmol-min.js'></script>" +
               "<style>html,body{margin:0;height:100%;background:#fff}#v{width:100%;height:100vh;position:relative}</style>" +
               "</head><body><div id='v'></div><script>" +
               "var changed=[" + resiList + "];" +
               "var viewer=$3Dmol.createViewer('v',{backgroundColor:'white'});" +
               "$3Dmol.download('pdb:6VXX',viewer,{},function(){" +
               "viewer.setStyle({},{cartoon:{color:'lightgray'}});" +
               "viewer.addStyle({resi:changed},{cartoon:{color:'red'}});" +
               "viewer.zoomTo();viewer.render();});" +
               "</script></body></html>";
    }

    private async void ShowMatrix_Click(object sender, RoutedEventArgs e)
    {
        bool pam = MatrixBox.SelectedIndex == 1;
        var m = pam ? Aligner.Pam250 : Aligner.Blosum62;
        const string letters = "ARNDCQEGHILKMFPSTWYV";
        var sb = new StringBuilder();
        sb.Append("   ");
        foreach (char c in letters) sb.Append(c.ToString().PadLeft(4));
        sb.Append('\n');
        foreach (char a in letters)
        {
            sb.Append(a).Append("  ");
            foreach (char b in letters) sb.Append(m[(a, b)].ToString().PadLeft(4));
            sb.Append('\n');
        }
        var dialog = new ContentDialog
        {
            Title = pam ? "PAM250 scoring matrix" : "BLOSUM62 scoring matrix",
            CloseButtonText = "Close",
            XamlRoot = XamlRoot,
            Content = new ScrollViewer
            {
                HorizontalScrollBarVisibility = ScrollBarVisibility.Auto,
                Content = new TextBlock { Text = sb.ToString(), FontFamily = new FontFamily("Consolas"), FontSize = 13 }
            }
        };
        await dialog.ShowAsync();
    }

    private async void LoadFile1_Click(object sender, RoutedEventArgs e)
    {
        var t = await PickFastaAsync();
        if (t != null) { Seq1Box.Text = t; RunAlign(); }
    }

    private async void LoadFile2_Click(object sender, RoutedEventArgs e)
    {
        var t = await PickFastaAsync();
        if (t != null) { Seq2Box.Text = t; RunAlign(); }
    }

    private async Task<string?> PickFastaAsync()
    {
        var picker = new FileOpenPicker();
        InitializeWithWindow.Initialize(picker, WindowNative.GetWindowHandle(App.MainWindow));
        picker.FileTypeFilter.Add(".fasta");
        picker.FileTypeFilter.Add(".fa");
        picker.FileTypeFilter.Add(".txt");
        picker.FileTypeFilter.Add("*");
        var file = await picker.PickSingleFileAsync();
        return file == null ? null : await FileIO.ReadTextAsync(file);
    }

    private void Copy_Click(object sender, RoutedEventArgs e)
    {
        var dp = new DataPackage();
        dp.SetText(_plainText);
        Clipboard.SetContent(dp);
    }

    private async void Save_Click(object sender, RoutedEventArgs e)
    {
        var picker = new FileSavePicker();
        InitializeWithWindow.Initialize(picker, WindowNative.GetWindowHandle(App.MainWindow));
        picker.FileTypeChoices.Add("Text file", new List<string> { ".txt" });
        picker.SuggestedFileName = "alignment";
        var file = await picker.PickSaveFileAsync();
        if (file != null) await FileIO.WriteTextAsync(file, _plainText);
    }

    // ---- run the alignment ----
    private async void RunAlign()
    {
        if (!_ready) return;
        int gen = ++_generation;

        string s1 = Clean(Seq1Box.Text);
        string s2 = Clean(Seq2Box.Text);
        if (s1.Length == 0 || s2.Length == 0)
        {
            StatusText.Text = "Type or load a sequence into both boxes.";
            AlignmentView.Blocks.Clear();
            return;
        }

        bool protein = TypeBox.SelectedIndex == 1;
        int mode = ModeBox.SelectedIndex; // 0 global, 1 local, 2 semi-global
        bool local = mode == 1;
        bool affine = GapModelBox.SelectedIndex == 1;
        Scorer sc = protein
            ? Aligner.FromMatrix(MatrixBox.SelectedIndex == 1 ? Aligner.Pam250 : Aligner.Blosum62)
            : Aligner.Dna(MatchSlider.Value, MismatchSlider.Value);
        double gap = GapSlider.Value;
        double open = OpenSlider.Value, ext = ExtendSlider.Value;

        long cells = (long)s1.Length * s2.Length;
        bool huge = cells > 9_000_000; // ~3000x3000; full tables would use too much memory
        StatusText.Text = huge ? "Long sequences: running a fast banded alignment..."
                               : (cells > 2_000_000 ? "Aligning..." : "");

        AlignResult? res = null;
        string note = "";
        await Task.Run(() =>
        {
            if (huge)
            {
                res = Aligner.Banded(s1, s2, sc, gap);
                if (mode != 0 || affine)
                    note = "These sequences are long, so a fast banded global alignment was used (local, semi-global, and affine are skipped at this size).";
                else
                    note = "Long sequences: aligned with the banded method.";
            }
            else if (affine && mode == 0) res = Aligner.GlobalAffine(s1, s2, sc, open, ext);
            else if (mode == 1) { res = Aligner.Local(s1, s2, sc, gap); if (affine) note = "Affine gaps apply to global alignment only, so linear was used here."; }
            else if (mode == 2) { res = Aligner.SemiGlobal(s1, s2, sc, gap); if (affine) note = "Affine gaps apply to global alignment only, so linear was used here."; }
            else res = Aligner.Global(s1, s2, sc, gap);
        });

        if (gen != _generation || res == null) return; // a newer request superseded this one

        _protein = protein;
        Render(res, local);
        StatusText.Text = note;
    }

    private void Render(AlignResult res, bool local)
    {
        ScoreText.Text = res.Score.ToString("0.##");
        IdentityText.Text = $"{Aligner.PercentIdentity(res.Row1, res.Row2):0.0}%";
        var (count, lengths) = Aligner.GapStats(res.Row1, res.Row2);
        GapsText.Text = count == 0 ? "0" : $"{count}  {Format(lengths)}";

        _plainText = BuildPlainText(res.Row1, res.Row2, local, res);
        BuildColoredView(res.Row1, res.Row2);
        DrawDotPlot(res.Row1.Replace("-", ""), res.Row2.Replace("-", ""));
        UpdateStats(res);
    }

    private void UpdateStats(AlignResult res)
    {
        string a = res.Row1.Replace("-", "");
        string b = res.Row2.Replace("-", "");
        int matches = 0, mism = 0;
        for (int i = 0; i < res.Row1.Length; i++)
        {
            if (res.Row1[i] == '-' || res.Row2[i] == '-') continue;
            if (res.Row1[i] == res.Row2[i]) matches++; else mism++;
        }
        string gc = _protein ? "" : $", GC {Gc(a):0}% / {Gc(b):0}%";
        StatsText.Text = $"Lengths {a.Length} and {b.Length}{gc}.  Aligned columns {res.Row1.Length}.  Matches {matches}, mismatches {mism}.";
    }

    private static double Gc(string s)
    {
        if (s.Length == 0) return 0;
        int gc = 0;
        foreach (char c in s) if (c == 'G' || c == 'C') gc++;
        return 100.0 * gc / s.Length;
    }

    private static string Format(List<int> lengths)
    {
        if (lengths.Count == 0) return "";
        var show = lengths.Count > 12 ? lengths.GetRange(0, 12) : lengths;
        return "[" + string.Join(",", show) + (lengths.Count > 12 ? ",..." : "") + "]";
    }

    private void BuildColoredView(string r1, string r2)
    {
        AlignmentView.Blocks.Clear();
        int shown = Math.Min(r1.Length, MaxViewColumns);
        const int width = 60;
        for (int start = 0; start < shown; start += width)
        {
            int end = Math.Min(start + width, shown);
            var p1 = new Paragraph();
            var pm = new Paragraph();
            var p2 = new Paragraph();
            for (int i = start; i < end; i++)
            {
                bool gap = r1[i] == '-' || r2[i] == '-';
                bool match = !gap && r1[i] == r2[i];
                var brush = gap ? GrayB : (match ? GreenB : RedB);
                char mc = gap ? ' ' : (match ? '|' : '.');
                p1.Inlines.Add(new Run { Text = r1[i].ToString(), Foreground = brush });
                pm.Inlines.Add(new Run { Text = mc.ToString(), Foreground = brush });
                p2.Inlines.Add(new Run { Text = r2[i].ToString(), Foreground = brush });
            }
            AlignmentView.Blocks.Add(p1);
            AlignmentView.Blocks.Add(pm);
            AlignmentView.Blocks.Add(p2);
            AlignmentView.Blocks.Add(new Paragraph());
        }
        if (r1.Length > MaxViewColumns)
        {
            var note = new Paragraph();
            note.Inlines.Add(new Run { Text = $"(showing first {MaxViewColumns} of {r1.Length} columns; use Save or Copy for the full alignment)", Foreground = GrayB });
            AlignmentView.Blocks.Add(note);
        }
    }

    private static string BuildPlainText(string r1, string r2, bool local, AlignResult res)
    {
        var sb = new StringBuilder();
        if (local)
            sb.AppendLine($"Local alignment. seq1[{res.Start1}:{res.End1}], seq2[{res.Start2}:{res.End2}]");
        const int width = 60;
        for (int start = 0; start < r1.Length; start += width)
        {
            int end = Math.Min(start + width, r1.Length);
            string a = r1.Substring(start, end - start);
            string b = r2.Substring(start, end - start);
            var mid = new StringBuilder();
            for (int i = 0; i < a.Length; i++)
                mid.Append(a[i] == '-' || b[i] == '-' ? ' ' : a[i] == b[i] ? '|' : '.');
            sb.AppendLine(a);
            sb.AppendLine(mid.ToString());
            sb.AppendLine(b);
            sb.AppendLine();
        }
        sb.AppendLine($"score = {res.Score:0.##}");
        sb.AppendLine($"percent identity = {Aligner.PercentIdentity(r1, r2):0.0}%");
        var (count, lengths) = Aligner.GapStats(r1, r2);
        sb.AppendLine($"gaps = {count}  lengths {Format(lengths)}");
        return sb.ToString();
    }

    private void DrawDotPlot(string s1, string s2)
    {
        DotCanvas.Children.Clear();
        double W = DotCanvas.Width, H = DotCanvas.Height;
        if (s1.Length == 0 || s2.Length == 0) return;
        if (s1.Length > 900 || s2.Length > 900)
        {
            DotCanvas.Children.Add(new TextBlock
            {
                Text = "Sequences are long, so the dot plot is skipped for speed.",
                Foreground = GrayB,
                TextWrapping = TextWrapping.Wrap,
                Width = W - 20
            });
            return;
        }
        int k = _protein ? 2 : 3;
        var index = new Dictionary<string, List<int>>();
        for (int j = 0; j + k <= s2.Length; j++)
        {
            string key = s2.Substring(j, k);
            if (!index.TryGetValue(key, out var list)) { list = new List<int>(); index[key] = list; }
            list.Add(j);
        }
        double sx = W / s2.Length, sy = H / s1.Length;
        for (int i = 0; i + k <= s1.Length; i++)
        {
            if (!index.TryGetValue(s1.Substring(i, k), out var js)) continue;
            foreach (int j in js)
            {
                var dot = new Rectangle { Width = 2, Height = 2, Fill = GreenB };
                Canvas.SetLeft(dot, j * sx);
                Canvas.SetTop(dot, i * sy);
                DotCanvas.Children.Add(dot);
            }
        }
    }
}
