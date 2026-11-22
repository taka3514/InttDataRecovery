#ifndef INTTWORDSHIFTCDBPRODUCTION_H
#define INTTWORDSHIFTCDBPRODUCTION_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>

class PHCompositeNode;
class TrkrClusterContainer;
class Gl1Packet;
class TFile;
class TTree;

class InttWordShiftCDBProduction : public SubsysReco
{
 public:
  explicit InttWordShiftCDBProduction(
      const std::string &name = "InttWordShiftCDBProduction");

  ~InttWordShiftCDBProduction() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;
  int Reset(PHCompositeNode *topNode) override;
  void Print(const std::string &what = "ALL") const override;

  void setOutFileName(const std::string &name);
  void setBinWidthSec(double sec);

 private:
  static constexpr int MAX_LAYER = 7;
  static constexpr int MAX_PHI = 16;
  static constexpr int MAX_SN = 2;

  static constexpr double BCO_SEC = 106e-9;

  struct BinData
  {
    uint64_t cluster_count = 0;
    double phi_size_sum = 0.0;

    void Reset()
    {
      cluster_count = 0;
      phi_size_sum = 0.0;
    }
  };

  void FinishBin(uint64_t end_bco);
  void ResetBin();

  std::string m_outfile_name = "intt_word_shift_stats.root";
  double m_bin_width_sec = 10.0;

  TFile *m_outfile = nullptr;
  TTree *m_tree_single = nullptr;
  TTree *m_tree_runrate = nullptr;

  TrkrClusterContainer *m_TrkrClusterContainer = nullptr;
  Gl1Packet *m_GL1Packet = nullptr;

  bool m_first_event = true;
  uint64_t m_bin_begin_bco = 0;
  uint64_t m_last_bco = 0;

  int m_segment = 0;

  BinData m_all[MAX_LAYER][MAX_PHI][MAX_SN];
  BinData m_tb0[MAX_LAYER][MAX_PHI][MAX_SN];

  // Single tree
  int m_layer = 0;
  int m_phi = 0;
  int m_sn = 0;

  int m_fs = 0;
  int m_fc = 0;

  uint64_t m_begin_bco = 0;
  uint64_t m_end_bco = 0;

  uint64_t m_cluster_count = 0;
  double m_cluster_rate = 0.0;
  double m_mean_phi_size = 0.0;
  double m_hit_rate = 0.0;

  uint64_t m_cluster_count_tb0 = 0;
  double m_cluster_rate_tb0 = 0.0;
  double m_mean_phi_size_tb0 = 0.0;
  double m_hit_rate_tb0 = 0.0;

  // RunRate tree
  int m_run = -1;

  uint64_t m_run_bco_count = 0;
  double m_run_time_sec = 0.0;

  uint64_t m_mbd_raw = 0;
  uint64_t m_mbd_live = 0;
  uint64_t m_mbd_scaled = 0;

  double m_mbd_scale_down = 0.0;
};

#endif