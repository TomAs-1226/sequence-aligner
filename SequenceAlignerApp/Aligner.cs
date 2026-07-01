using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SequenceAlignerApp;

/// <summary>Scores a pair of aligned symbols.</summary>
public delegate double Scorer(char a, char b);

public sealed record AlignResult(
    string Row1, string Row2, double Score,
    int Start1 = 0, int End1 = 0, int Start2 = 0, int End2 = 0);

/// <summary>
/// C# port of the Python alignment engine: global and local alignment with
/// linear or affine gaps, DNA or protein scoring, plus readouts and translation.
/// Same algorithms the notebook uses.
/// </summary>
public static class Aligner
{
    private const double Eps = 1e-9;

    // ---- scorers ----
    public static Scorer Dna(double match, double mismatch) =>
        (a, b) => a == b ? match : -mismatch;

    public static Scorer FromMatrix(Dictionary<(char, char), int> m) => (a, b) =>
    {
        if (m.TryGetValue((a, b), out var v)) return v;
        char aa = m.ContainsKey((a, 'A')) ? a : 'X';
        char bb = m.ContainsKey(('A', b)) ? b : 'X';
        return m[(aa, bb)];
    };

    // ---- global alignment (linear gap) ----
    public static AlignResult Global(string s1, string s2, Scorer sc, double gap)
    {
        int n = s1.Length, m = s2.Length;
        var S = new double[n + 1, m + 1];
        for (int i = 1; i <= n; i++) S[i, 0] = S[i - 1, 0] - gap;
        for (int j = 1; j <= m; j++) S[0, j] = S[0, j - 1] - gap;
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
            {
                double diag = S[i - 1, j - 1] + sc(s1[i - 1], s2[j - 1]);
                double up = S[i - 1, j] - gap;
                double left = S[i, j - 1] - gap;
                S[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        var (r1, r2, _, _) = Traceback(s1, s2, sc, gap, S, n, m, local: false);
        return new AlignResult(r1, r2, S[n, m]);
    }

    // ---- local alignment (linear gap) ----
    public static AlignResult Local(string s1, string s2, Scorer sc, double gap)
    {
        int n = s1.Length, m = s2.Length;
        var S = new double[n + 1, m + 1];
        double best = 0; int bi = 0, bj = 0;
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
            {
                double diag = S[i - 1, j - 1] + sc(s1[i - 1], s2[j - 1]);
                double up = S[i - 1, j] - gap;
                double left = S[i, j - 1] - gap;
                double v = Math.Max(0, Math.Max(diag, Math.Max(up, left)));
                S[i, j] = v;
                if (v > best) { best = v; bi = i; bj = j; }
            }
        var (r1, r2, si, sj) = Traceback(s1, s2, sc, gap, S, bi, bj, local: true);
        return new AlignResult(r1, r2, best, si, bi, sj, bj);
    }

    private static (string, string, int, int) Traceback(
        string s1, string s2, Scorer sc, double gap, double[,] S, int i, int j, bool local)
    {
        var a = new StringBuilder();
        var b = new StringBuilder();
        while ((i > 0 || j > 0) && !(local && S[i, j] <= Eps))
        {
            if (i > 0 && j > 0 && Math.Abs(S[i, j] - (S[i - 1, j - 1] + sc(s1[i - 1], s2[j - 1]))) < Eps)
            { a.Append(s1[i - 1]); b.Append(s2[j - 1]); i--; j--; }
            else if (i > 0 && Math.Abs(S[i, j] - (S[i - 1, j] - gap)) < Eps)
            { a.Append(s1[i - 1]); b.Append('-'); i--; }
            else { a.Append('-'); b.Append(s2[j - 1]); j--; }
        }
        return (Rev(a), Rev(b), i, j);
    }

    // ---- global alignment (affine gap, Gotoh) ----
    public static AlignResult GlobalAffine(string s1, string s2, Scorer sc, double open, double extend)
    {
        int n = s1.Length, m = s2.Length;
        double NEG = double.NegativeInfinity;
        var M = new double[n + 1, m + 1];
        var Ix = new double[n + 1, m + 1];
        var Iy = new double[n + 1, m + 1];
        for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) { M[i, j] = NEG; Ix[i, j] = NEG; Iy[i, j] = NEG; }
        M[0, 0] = 0;
        for (int i = 1; i <= n; i++) Ix[i, 0] = -open - (i - 1) * extend;
        for (int j = 1; j <= m; j++) Iy[0, j] = -open - (j - 1) * extend;
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
            {
                double sub = sc(s1[i - 1], s2[j - 1]);
                M[i, j] = Math.Max(M[i - 1, j - 1], Math.Max(Ix[i - 1, j - 1], Iy[i - 1, j - 1])) + sub;
                Ix[i, j] = Math.Max(M[i - 1, j] - open, Ix[i - 1, j] - extend);
                Iy[i, j] = Math.Max(M[i, j - 1] - open, Iy[i, j - 1] - extend);
            }
        // traceback
        var a = new StringBuilder(); var b = new StringBuilder();
        int ii = n, jj = m;
        string state = "M";
        double best = Math.Max(M[n, m], Math.Max(Ix[n, m], Iy[n, m]));
        if (best == Ix[n, m]) state = "Ix"; else if (best == Iy[n, m]) state = "Iy";
        while (ii > 0 || jj > 0)
        {
            if (state == "M")
            {
                a.Append(s1[ii - 1]); b.Append(s2[jj - 1]);
                double sub = sc(s1[ii - 1], s2[jj - 1]);
                double prev = M[ii, jj] - sub;
                ii--; jj--;
                if (Math.Abs(prev - M[ii, jj]) < Eps) state = "M";
                else if (Math.Abs(prev - Ix[ii, jj]) < Eps) state = "Ix";
                else state = "Iy";
            }
            else if (state == "Ix")
            {
                a.Append(s1[ii - 1]); b.Append('-');
                state = (ii >= 1 && Math.Abs(Ix[ii, jj] - (M[ii - 1, jj] - open)) < Eps) ? "M" : "Ix";
                ii--;
            }
            else
            {
                a.Append('-'); b.Append(s2[jj - 1]);
                state = (jj >= 1 && Math.Abs(Iy[ii, jj] - (M[ii, jj - 1] - open)) < Eps) ? "M" : "Iy";
                jj--;
            }
        }
        return new AlignResult(Rev(a), Rev(b), best);
    }

    // ---- readouts ----
    public static double PercentIdentity(string r1, string r2)
    {
        int matches = 0, aligned = 0;
        for (int i = 0; i < r1.Length; i++)
        {
            if (r1[i] == '-' || r2[i] == '-') continue;
            aligned++;
            if (r1[i] == r2[i]) matches++;
        }
        return aligned == 0 ? 0 : 100.0 * matches / aligned;
    }

    public static (int count, List<int> lengths) GapStats(string r1, string r2)
    {
        var lengths = new List<int>();
        int run = 0;
        for (int i = 0; i < r1.Length; i++)
        {
            if (r1[i] == '-' || r2[i] == '-') run++;
            else { if (run > 0) lengths.Add(run); run = 0; }
        }
        if (run > 0) lengths.Add(run);
        return (lengths.Count, lengths);
    }

    public static double AlignmentScore(string r1, string r2, Scorer sc, double gap)
    {
        double total = 0;
        for (int i = 0; i < r1.Length; i++)
        {
            if (r1[i] == '-' || r2[i] == '-') total -= gap;
            else total += sc(r1[i], r2[i]);
        }
        return total;
    }

    // ---- translation ----
    public static string Translate(string dna)
    {
        dna = dna.ToUpperInvariant().Replace("U", "T");
        var sb = new StringBuilder();
        for (int i = 0; i + 3 <= dna.Length; i += 3)
        {
            string codon = dna.Substring(i, 3);
            if (!Codons.TryGetValue(codon, out var aa)) aa = 'X';
            if (aa == '*') break;
            sb.Append(aa);
        }
        return sb.ToString();
    }

    private static string Rev(StringBuilder sb)
    {
        var arr = sb.ToString().ToCharArray();
        Array.Reverse(arr);
        return new string(arr);
    }

    // ---- data: codon table + matrices ----
    private static readonly Dictionary<string, char> Codons = BuildCodons();

    private static Dictionary<string, char> BuildCodons()
    {
        string bases = "TCAG";
        // amino acids in the standard TCAG x TCAG x TCAG order
        string aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        var d = new Dictionary<string, char>();
        int k = 0;
        foreach (char b1 in bases)
            foreach (char b2 in bases)
                foreach (char b3 in bases)
                    d[$"{b1}{b2}{b3}"] = aas[k++];
        return d;
    }

    public static readonly Dictionary<(char, char), int> Blosum62 = ParseMatrix(Blosum62Text);
    public static readonly Dictionary<(char, char), int> Pam250 = ParseMatrix(Pam250Text);

    private static Dictionary<(char, char), int> ParseMatrix(string text)
    {
        var d = new Dictionary<(char, char), int>();
        var lines = text.Trim().Replace("\r", "").Split('\n');
        var cols = lines[0].Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
        for (int r = 1; r < lines.Length; r++)
        {
            var parts = lines[r].Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
            char row = parts[0][0];
            for (int c = 0; c < cols.Length; c++)
                d[(row, cols[c][0])] = int.Parse(parts[c + 1]);
        }
        return d;
    }

    private const string Blosum62Text = @"
   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
A   4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0  -2  -1   0  -4
R  -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3  -1   0  -1  -4
N  -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3   3   0  -1  -4
D  -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3   4   1  -1  -4
C   0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1  -3  -3  -2  -4
Q  -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2   0   3  -1  -4
E  -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
G   0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3  -1  -2  -1  -4
H  -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3   0   0  -1  -4
I  -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3  -3  -3  -1  -4
L  -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1  -4  -3  -1  -4
K  -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2   0   1  -1  -4
M  -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1  -3  -1  -1  -4
F  -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1  -3  -3  -1  -4
P  -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2  -2  -1  -2  -4
S   1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -2   0   0   0  -4
T   0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0  -1  -1   0  -4
W  -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3  -4  -3  -2  -4
Y  -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1  -3  -2  -1  -4
V   0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -2   0  -3  -1   4  -3  -2  -1  -4
B  -2  -1   3   4  -3   0   1  -1   0  -3  -4   0  -3  -3  -2   0  -1  -4  -3  -3   4   1  -1  -4
Z  -1   0   0   1  -3   3   4  -2   0  -3  -3   1  -1  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
X   0  -1  -1  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2   0   0  -2  -1  -1  -1  -1  -1  -4
*  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4   1";

    private const string Pam250Text = @"
   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
A   2  -2   0   0  -2   0   0   1  -1  -1  -2  -1  -1  -3   1   1   1  -6  -3   0   0   0   0  -8
R  -2   6   0  -1  -4   1  -1  -3   2  -2  -3   3   0  -4   0   0  -1   2  -4  -2  -1   0  -1  -8
N   0   0   2   2  -4   1   1   0   2  -2  -3   1  -2  -3   0   1   0  -4  -2  -2   2   1   0  -8
D   0  -1   2   4  -5   2   3   1   1  -2  -4   0  -3  -6  -1   0   0  -7  -4  -2   3   3  -1  -8
C  -2  -4  -4  -5  12  -5  -5  -3  -3  -2  -6  -5  -5  -4  -3   0  -2  -8   0  -2  -4  -5  -3  -8
Q   0   1   1   2  -5   4   2  -1   3  -2  -2   1  -1  -5   0  -1  -1  -5  -4  -2   1   3  -1  -8
E   0  -1   1   3  -5   2   4   0   1  -2  -3   0  -2  -5  -1   0   0  -7  -4  -2   3   3  -1  -8
G   1  -3   0   1  -3  -1   0   5  -2  -3  -4  -2  -3  -5   0   1   0  -7  -5  -1   0   0  -1  -8
H  -1   2   2   1  -3   3   1  -2   6  -2  -2   0  -2  -2   0  -1  -1  -3   0  -2   1   2  -1  -8
I  -1  -2  -2  -2  -2  -2  -2  -3  -2   5   2  -2   2   1  -2  -1   0  -5  -1   4  -2  -2  -1  -8
L  -2  -3  -3  -4  -6  -2  -3  -4  -2   2   6  -3   4   2  -3  -3  -2  -2  -1   2  -3  -3  -1  -8
K  -1   3   1   0  -5   1   0  -2   0  -2  -3   5   0  -5  -1   0   0  -3  -4  -2   1   0  -1  -8
M  -1   0  -2  -3  -5  -1  -2  -3  -2   2   4   0   6   0  -2  -2  -1  -4  -2   2  -2  -2  -1  -8
F  -3  -4  -3  -6  -4  -5  -5  -5  -2   1   2  -5   0   9  -5  -3  -3   0   7  -1  -4  -5  -2  -8
P   1   0   0  -1  -3   0  -1   0   0  -2  -3  -1  -2  -5   6   1   0  -6  -5  -1  -1   0  -1  -8
S   1   0   1   0   0  -1   0   1  -1  -1  -3   0  -2  -3   1   2   1  -2  -3  -1   0   0   0  -8
T   1  -1   0   0  -2  -1   0   0  -1   0  -2   0  -1  -3   0   1   3  -5  -3   0   0  -1   0  -8
W  -6   2  -4  -7  -8  -5  -7  -7  -3  -5  -2  -3  -4   0  -6  -2  -5  17   0  -6  -5  -6  -4  -8
Y  -3  -4  -2  -4   0  -4  -4  -5   0  -1  -1  -4  -2   7  -5  -3  -3   0  10  -2  -3  -4  -2  -8
V   0  -2  -2  -2  -2  -2  -2  -1  -2   4   2  -2   2  -1  -1  -1   0  -6  -2   4  -2  -2  -1  -8
B   0  -1   2   3  -4   1   3   0   1  -2  -3   1  -2  -4  -1   0   0  -5  -3  -2   3   2  -1  -8
Z   0   0   1   3  -5   3   3   0   2  -2  -3   0  -2  -5   0   0  -1  -6  -4  -2   2   3  -1  -8
X   0  -1   0  -1  -3  -1  -1  -1  -1  -1  -1  -1  -1  -2  -1   0   0  -4  -2  -1  -1  -1  -1  -8
*  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8   1";
}
