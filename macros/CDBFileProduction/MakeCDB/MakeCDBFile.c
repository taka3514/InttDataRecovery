#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TLine.h>
#include <TDirectory.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLatex.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr double BCO_SEC = 106e-9;

struct Key {
  int layer, phi, sn;

  bool operator<(const Key& b) const {
    return std::tie(layer, phi, sn) <
           std::tie(b.layer, b.phi, b.sn);
  }
};

struct Bin {
  ULong64_t begin, end, count;
  double phi, rate;
  int fs, fc;
};

struct Range {
  ULong64_t begin, end;
};

double median(std::vector<double> v) {
  if (v.empty())
    return std::numeric_limits<double>::quiet_NaN();

  std::sort(v.begin(), v.end());

  size_t n = v.size();

  return n % 2
    ? v[n / 2]
    : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

}

// Usage:
// root -l -b -q 'CalibFileProductionFromStats.C(79516,0,40000,1.55,0.25,1.0)'

void MakeCDBFile(
  int run = 81000,
  int seg_min = 0,
  int seg_max = 40000,
  double phi_size_cut = 1.5,
  double center_tolerance = 0.2,
  double connect_gap_sec = 1.0,
  const char* output_file = "")
{
  if (seg_min < 0 ||
      seg_max < seg_min ||
      phi_size_cut <= 0 ||
      center_tolerance < 0 ||
      center_tolerance >= 1 ||
      connect_gap_sec < 0) {

    std::cerr << "ERROR: invalid arguments"
              << std::endl;
    return;
  }

  const std::string dir = Form(
    "/sphenix/tg/tg01/commissioning/INTT/work/takakiku/"
    "GitHub_OutPut/ScanTrkrCluster/%d",
    run
  );

  // const std::string dir = Form(
  //   "/sphenix/u/takakiku/user_takakiku/GitHub/"
  //   "InttDataRecovery/macros/CDBFileProduction/MakeCDB"
  // );

  // const std::string out =
  //   (output_file && *output_file)
  //   ? output_file
  //   : Form("%s/calib_%08d.root", dir.c_str(), run);

  const std::string cur_dir =
    "/sphenix/tg/tg01/commissioning/INTT/work/takakiku/"
    "GitHub_OutPut/WordShiftCDB";

  const std::string out =
    (output_file && *output_file)
    ? output_file
    : Form("%s/calib_%08d.root", cur_dir.c_str(), run);

  std::map<Key, std::vector<Bin>> bins;
  std::map<Key, std::pair<int, int>> mapping;

  ULong64_t first_bco =
    std::numeric_limits<ULong64_t>::max();

  ULong64_t last_bco = 0;

  int nfiles = 0;
  int invalid = 0;

  for (int seg = seg_min; seg <= seg_max; ++seg) {

    const std::string filename = Form(
      "%s/stats_%08d_%05d.root",
      dir.c_str(),
      run,
      seg
    );

    if (gSystem->AccessPathName(filename.c_str()))
      continue;

    std::unique_ptr<TFile> f(
      TFile::Open(filename.c_str(), "READ")
    );

    if (!f || f->IsZombie()) {
      ++invalid;
      continue;
    }

    TTree* t = nullptr;
    f->GetObject("Single", t);

    if (!t) {
      ++invalid;
      std::cerr << "Missing Single: "
                << filename << std::endl;
      continue;
    }

    const char* required[] = {
      "layer",
      "phi",
      "sn",
      "felix_server",
      "felix_channel",
      "begin_bco",
      "end_bco",
      "cluster_count",
      "mean_phi_size",
      "cluster_rate_tb0",
      "hit_rate_tb0"
    };

    bool valid = true;

    for (auto name : required) {
      if (!t->GetBranch(name)) {
        std::cerr << "Missing "
                  << name
                  << " in "
                  << filename
                  << std::endl;

        valid = false;
      }
    }

    if (!valid) {
      ++invalid;
      continue;
    }

    int layer = 0;
    int phi = 0;
    int sn = 0;
    int fs = 0;
    int fc = 0;

    ULong64_t begin = 0;
    ULong64_t end = 0;
    ULong64_t count = 0;

    double mean = 0;
    double rate = 0;
    double hit_rate = 0;

    t->SetBranchAddress("layer", &layer);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("sn", &sn);

    t->SetBranchAddress("felix_server", &fs);
    t->SetBranchAddress("felix_channel", &fc);

    t->SetBranchAddress("begin_bco", &begin);
    t->SetBranchAddress("end_bco", &end);

    t->SetBranchAddress("cluster_count", &count);
    t->SetBranchAddress("mean_phi_size", &mean);

    t->SetBranchAddress("cluster_rate_tb0", &rate);
    t->SetBranchAddress("hit_rate_tb0", &hit_rate);

    for (Long64_t i = 0; i < t->GetEntries(); ++i) {

      t->GetEntry(i);

      if (layer < 3 ||
          layer > 6 ||
          phi < 0 ||
          phi >= ((layer < 5) ? 12 : 16) ||
          sn < 0 ||
          sn > 1 ||
          end <= begin)
        continue;

      Key k{layer, phi, sn};

      auto it = mapping.find(k);

      if (it != mapping.end() &&
          it->second != std::make_pair(fs, fc)) {

        std::cerr
          << "ERROR: inconsistent FELIX mapping in "
          << filename
          << std::endl;

        return;
      }

      mapping[k] = {fs, fc};

      bins[k].push_back({
        begin,
        end,
        count,
        mean,
        rate,
        fs,
        fc
      });

      first_bco = std::min(first_bco, begin);
      last_bco = std::max(last_bco, end);
    }

    ++nfiles;

    if (nfiles % 1000 == 0) {
      std::cout
        << "Processed "
        << nfiles
        << " files (current segment = "
        << seg
        << ")"
        << std::endl;
    }
  }

  if (!nfiles ||
      first_bco ==
        std::numeric_limits<ULong64_t>::max()) {

    std::cerr
      << "ERROR: no valid stats input in "
      << dir
      << std::endl;

    return;
  }

  std::unique_ptr<TFile> fout(
    TFile::Open(out.c_str(), "RECREATE")
  );

  if (!fout || fout->IsZombie()) {
    std::cerr
      << "ERROR: cannot create "
      << out
      << std::endl;

    return;
  }

  TTree tree(
    "Multiple",
    "Bad periods and half-ladder classification from stats"
  );

  int o_fs = -1;
  int o_fc = -1;
  int o_layer = -1;
  int o_phi = -1;
  int o_sn = -1;

  std::vector<ULong64_t> o_bad_begin;
  std::vector<ULong64_t> o_bad_end;
  std::vector<int> o_bad_type;

  double o_normal_center = 0;
  double o_bad_center = 0;
  double o_hit_rate_ratio = 0;

  tree.Branch(
    "felix_server",
    &o_fs,
    "felix_server/I"
  );

  tree.Branch(
    "felix_channel",
    &o_fc,
    "felix_channel/I"
  );

  tree.Branch("layer", &o_layer, "layer/I");
  tree.Branch("phi", &o_phi, "phi/I");
  tree.Branch("sn", &o_sn, "sn/I");

  tree.Branch("bad_begin", &o_bad_begin);
  tree.Branch("bad_end", &o_bad_end);
  tree.Branch("bad_type", &o_bad_type);

  tree.Branch(
    "normal_center",
    &o_normal_center,
    "normal_center/D"
  );

  tree.Branch(
    "bad_center",
    &o_bad_center,
    "bad_center/D"
  );

  tree.Branch(
    "hit_rate_ratio",
    &o_hit_rate_ratio,
    "hit_rate_ratio/D"
  );

  TDirectory* plots = fout->mkdir("Plots");

  const double duration =
    (last_bco - first_bco) * BCO_SEC;

  long long nbad = 0;

  for (int layer = 3; layer <= 6; ++layer) {
    for (int phi = 0;
         phi < ((layer < 5) ? 12 : 16);
         ++phi) {
      for (int sn = 0; sn < 2; ++sn) {

        Key k{layer, phi, sn};

        auto& v = bins[k];

        std::sort(
          v.begin(),
          v.end(),
          [](const Bin& a, const Bin& b) {
            return std::tie(a.begin, a.end) <
                   std::tie(b.begin, b.end);
          }
        );

        o_layer = layer;
        o_phi = phi;
        o_sn = sn;
        o_fs = -1;
        o_fc = -1;

        if (mapping.count(k)) {
          o_fs = mapping[k].first;
          o_fc = mapping[k].second;
        }

        o_bad_begin.clear();
        o_bad_end.clear();
        o_bad_type.clear();

        std::vector<double> good_rates;
        std::vector<double> bad_rates;
        std::vector<Range> periods;

        bool active = false;

        ULong64_t bad_start = 0;
        ULong64_t bad_end = 0;

        // A gap with no input bins is NOT
        // automatically classified as bad.
        // Adjacent bad bins may be connected
        // only within connect_gap_sec.

        const ULong64_t gap_bco =
          static_cast<ULong64_t>(
            std::llround(connect_gap_sec / BCO_SEC)
          );

        for (const auto& b : v) {

          if (b.count == 0 ||
              !std::isfinite(b.phi) ||
              b.phi <= 0)
            continue;

          if (b.rate > 0 && std::isfinite(b.rate)) {

            if (b.phi < phi_size_cut)
              bad_rates.push_back(b.rate);
            else
              good_rates.push_back(b.rate);
          }

          if (b.phi < phi_size_cut) {

            if (!active) {
              bad_start = b.begin;
              bad_end = b.end;
              active = true;
            }
            else if (
              b.begin <= bad_end ||
              b.begin - bad_end <= gap_bco
            ) {
              bad_end = std::max(bad_end, b.end);
            }
            else {
              periods.push_back({
                bad_start,
                bad_end
              });

              bad_start = b.begin;
              bad_end = b.end;
            }
          }
          else if (active) {

            periods.push_back({
              bad_start,
              bad_end
            });

            active = false;
          }
        }

        if (active) {
          periods.push_back({
            bad_start,
            bad_end
          });
        }

        o_normal_center = median(good_rates);
        o_bad_center = median(bad_rates);

        o_hit_rate_ratio =
          std::numeric_limits<double>::quiet_NaN();

        int type = -1;

        if (std::isfinite(o_normal_center) &&
            o_normal_center > 0 &&
            std::isfinite(o_bad_center)) {

          o_hit_rate_ratio =
            o_bad_center / o_normal_center;

          if (o_bad_center > o_normal_center)
            type = 0;  // Enhancement
          else if (o_bad_center < o_normal_center)
            type = 1;  // Dropping
          else
            type = -1; // Equal
        }

        for (const auto& p : periods) {

          o_bad_begin.push_back(p.begin);
          o_bad_end.push_back(p.end);
          o_bad_type.push_back(type);

          ++nbad;
        }

        tree.Fill();

        if (v.empty())
          continue;

        plots->cd();

        const std::string stem = Form(
          "layer%d_phi%d_sn%d",
          layer,
          phi,
          sn
        );

        const int nb = std::max(
          1,
          static_cast<int>(v.size())
        );

        // Variable bin edges preserve original
        // time-bin widths, including gaps.
        // For irregular/overlapping segment bins,
        // use a fine fixed grid instead.

        std::vector<double> edges;

        bool regular = true;

        for (size_t i = 0; i < v.size(); ++i) {

          double a =
            (v[i].begin - first_bco) * BCO_SEC;

          double b =
            (v[i].end - first_bco) * BCO_SEC;

          if (i && a < edges.back() - 1e-8)
            regular = false;

          if (i && a > edges.back() + 1e-8)
            regular = false;

          if (i == 0)
            edges.push_back(a);

          edges.push_back(b);
        }

        std::unique_ptr<TH1D> hphi;
        std::unique_ptr<TH1D> hrate;

        if (regular) {

          hphi.reset(
            new TH1D(
              ("h_phi_" + stem).c_str(),
              ";Time since first BCO [sec];Mean cluster phi size (TB=any)",
              nb,
              edges.data()
            )
          );

          hrate.reset(
            new TH1D(
              ("h_rate_" + stem).c_str(),
              ";Time since first BCO [sec];Cluster yield rate (TB=0) [Hz]",
              nb,
              edges.data()
            )
          );

          for (size_t i = 0; i < v.size(); ++i) {

            hphi->SetBinContent(
              i + 1,
              v[i].count ? v[i].phi : 0
            );

            hrate->SetBinContent(
              i + 1,
              v[i].rate
            );
          }
        }
        else {

          // TH1D is a binned display:
          // bin centers represent sampled times.

         // Fixed-width bins: approximately 1 second per bin
        int ngrid = std::max(
            1,
            static_cast<int>(std::ceil(duration / 1.0))
        );

          hphi.reset(
            new TH1D(
              ("h_phi_" + stem).c_str(),
              ";Time since first BCO [sec];Mean cluster phi size (TB=any)",
              ngrid,
              0,
              duration
            )
          );

          hrate.reset(
            new TH1D(
              ("h_rate_" + stem).c_str(),
              ";Time since first BCO [sec];Cluster yield rate (TB=0) [Hz]",
              ngrid,
              0,
              duration
            )
          );

          for (const auto& b : v) {

            double x =
              ((b.begin - first_bco) +
               (b.end - b.begin) / 2.0) * BCO_SEC;

            int ib = hphi->FindBin(x);

            if (ib >= 1 && ib <= ngrid) {

              hphi->SetBinContent(
                ib,
                b.count ? b.phi : 0
              );

              hrate->SetBinContent(
                ib,
                b.rate
              );
            }
          }
        }

        hphi->SetDirectory(nullptr);
        hrate->SetDirectory(nullptr);

        TCanvas c(
          ("c_" + stem).c_str(),
          stem.c_str(),
          1300,
          850
        );

        c.Divide(1, 2);

        for (int pad = 1; pad <= 2; ++pad) {

          c.cd(pad);

          // ==================================
          // Disable background grid lines
          // ==================================
          gPad->SetGridx(0);
          gPad->SetGridy(0);

          TH1D* h =
            pad == 1 ? hphi.get() : hrate.get();

          h->SetLineColor(kTeal + 2);
          h->SetLineWidth(2);

          h->SetMinimum(0);

          double ymax = std::max(
            1.0,
            h->GetMaximum() * 1.2
          );

          if (pad == 1)
            ymax = std::max(
              ymax,
              phi_size_cut * 1.25
            );

          // Keep the cluster-yield median lines
          // within the displayed y range.

          if (pad == 2 &&
              std::isfinite(o_normal_center) &&
              o_normal_center > 0) {

            ymax = std::max(
              ymax,
              o_normal_center * 1.2
            );
          }

          if (pad == 2 &&
              std::isfinite(o_bad_center) &&
              o_bad_center >= 0) {

            ymax = std::max(
              ymax,
              o_bad_center * 1.2
            );
          }

          h->SetMaximum(ymax);
          h->Draw("HIST");

          // ==================================
          // Mean cluster phi size threshold
          // Black SOLID line
          // ==================================

          if (pad == 1) {

            TLine cut(
              0,
              phi_size_cut,
              duration,
              phi_size_cut
            );

            cut.SetLineColor(kBlack);
            cut.SetLineWidth(2);
            cut.SetLineStyle(1);
            cut.DrawClone();
          }

          // ==================================
          // Cluster-yield reference levels
          // Normal / Failure medians
          // ==================================

          if (pad == 2 &&
              std::isfinite(o_normal_center) &&
              o_normal_center > 0) {

            // Normal median:
            // Black SOLID line

            TLine normal(
              0,
              o_normal_center,
              duration,
              o_normal_center
            );

            normal.SetLineColor(kBlack);
            normal.SetLineStyle(1);
            normal.SetLineWidth(2);
            normal.DrawClone();

            // Failure median:
            // Magenta solid line

            if (std::isfinite(o_bad_center)) {

              TLine failure(
                0,
                o_bad_center,
                duration,
                o_bad_center
              );

              failure.SetLineColor(kMagenta + 2);
              failure.SetLineStyle(1);
              failure.SetLineWidth(3);
              failure.DrawClone();
            }
          }

          // ==================================
          // Bad-period boundaries
          // Red = start
          // Blue = end
          // ==================================

          for (const auto& p : periods) {

            double a =
              (p.begin - first_bco) * BCO_SEC;

            double b =
              (p.end - first_bco) * BCO_SEC;

            TLine start(a, 0, a, ymax);

            start.SetLineColor(kRed);
            start.SetLineWidth(2);
            start.DrawClone();

            TLine stop(b, 0, b, ymax);

            stop.SetLineColor(kBlue);
            stop.SetLineWidth(2);
            stop.DrawClone();
          }
        }

        // Show the existing bad_type classification
        // and its yield ratio on the canvas.

        const char* classification =
          periods.empty()
            ? "NORMAL"
            : (
                type == 0
                  ? "ENHANCEMENT"
                  : (
                      type == 1
                        ? "DROPPING"
                        : "EQUAL / UNKNOWN"
                    )
              );

        const int classification_color =
          periods.empty()
            ? kGreen + 2
            : (
                type == 0
                  ? kRed + 1
                  : (
                      type == 1
                        ? kBlue + 1
                        : kGray + 2
                    )
              );

        c.cd(1);

        TLatex label;

        label.SetNDC();
        label.SetTextFont(42);
        label.SetTextSize(0.040);
        label.SetTextColor(classification_color);

        if (std::isfinite(o_hit_rate_ratio)) {

          label.DrawLatex(
            0.15,
            0.85,
            Form(
              "%s  (Cluster yield ratio = %.3f)",
              classification,
              o_hit_rate_ratio
            )
          );
        }
        else {

          label.DrawLatex(
            0.15,
            0.85,
            classification
          );
        }

        c.cd(2);

        TLatex yield_label;

        yield_label.SetNDC();
        yield_label.SetTextFont(42);
        yield_label.SetTextSize(0.030);

        yield_label.SetTextColor(kBlack);

        if (std::isfinite(o_normal_center)) {

          yield_label.DrawLatex(
            0.13,
            0.85,
            Form(
              "Normal median: %.3g",
              o_normal_center
            )
          );
        }

        yield_label.SetTextColor(kMagenta + 2);

        if (std::isfinite(o_bad_center)) {

          yield_label.DrawLatex(
            0.13,
            0.79,
            Form(
              "Failure median: %.3g",
              o_bad_center
            )
          );
        }

        plots->cd();

        hphi->Write();
        hrate->Write();
        c.Write();
      }
    }
  }

  fout->cd();

  tree.Write();
  fout->Write();
  fout->Close();

  std::cout
    << "Saved " << out
    << " | files=" << nfiles
    << " invalid=" << invalid
    << " | half-ladders=" << tree.GetEntries()
    << " | bad periods=" << nbad
    << " | first BCO=" << first_bco
    << std::endl;
}