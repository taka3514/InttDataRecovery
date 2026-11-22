#include "InttWordShiftCDBProduction.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/InttDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <ffarawobjects/Gl1Packet.h>
#include <calotrigger/TriggerRunInfo.h>

#include <intt/InttMapping.h>

#include <TFile.h>
#include <TTree.h>

#include <iostream>

//____________________________________________________________________________..
InttWordShiftCDBProduction::InttWordShiftCDBProduction(
    const std::string &name)
  : SubsysReco(name)
{
  std::cout
      << "InttWordShiftCDBProduction::InttWordShiftCDBProduction()"
      << std::endl;
}

//____________________________________________________________________________..
InttWordShiftCDBProduction::~InttWordShiftCDBProduction()
{
  std::cout
      << "InttWordShiftCDBProduction::~InttWordShiftCDBProduction()"
      << std::endl;
}

//____________________________________________________________________________..
int InttWordShiftCDBProduction::Init(
    [[maybe_unused]] PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttWordShiftCDBProduction::InitRun(
    PHCompositeNode *topNode)
{
  // ============================================================
  // Input nodes
  // ============================================================

  m_TrkrClusterContainer =
      findNode::getClass<TrkrClusterContainer>(
          topNode, "TRKR_CLUSTER");

  if (!m_TrkrClusterContainer)
  {
    std::cerr << "ERROR: TRKR_CLUSTER not found" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_GL1Packet =
      findNode::getClass<Gl1Packet>(
          topNode, "GL1RAWHIT");

  if (!m_GL1Packet)
  {
    std::cerr << "ERROR: GL1RAWHIT not found" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_bin_width_sec <= 0.0)
  {
    std::cerr << "ERROR: bin width must be > 0" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // ============================================================
  // Output ROOT file
  // ============================================================

  m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");

  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << "ERROR: cannot create output file" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // ============================================================
  // Single tree
  // One entry = one time bin x one half-ladder
  // ============================================================

  m_tree_single = new TTree(
      "Single",
      "INTT half-ladder time-bin statistics");

  m_tree_single->Branch("segment", &m_segment, "segment/I");

  m_tree_single->Branch("begin_bco", &m_begin_bco, "begin_bco/l");
  m_tree_single->Branch("end_bco", &m_end_bco, "end_bco/l");

  m_tree_single->Branch("layer", &m_layer, "layer/I");
  m_tree_single->Branch("phi", &m_phi, "phi/I");
  m_tree_single->Branch("sn", &m_sn, "sn/I");

  m_tree_single->Branch("felix_server", &m_fs, "felix_server/I");
  m_tree_single->Branch("felix_channel", &m_fc, "felix_channel/I");

  // TB = any
  m_tree_single->Branch(
      "cluster_count", &m_cluster_count, "cluster_count/l");
  m_tree_single->Branch(
      "cluster_rate", &m_cluster_rate, "cluster_rate/D");
  m_tree_single->Branch(
      "mean_phi_size", &m_mean_phi_size, "mean_phi_size/D");
  m_tree_single->Branch(
      "hit_rate", &m_hit_rate, "hit_rate/D");

  // TB = 0
  m_tree_single->Branch(
      "cluster_count_tb0", &m_cluster_count_tb0, "cluster_count_tb0/l");
  m_tree_single->Branch(
      "cluster_rate_tb0", &m_cluster_rate_tb0, "cluster_rate_tb0/D");
  m_tree_single->Branch(
      "mean_phi_size_tb0", &m_mean_phi_size_tb0, "mean_phi_size_tb0/D");
  m_tree_single->Branch(
      "hit_rate_tb0", &m_hit_rate_tb0, "hit_rate_tb0/D");

  // ============================================================
  // RunRate tree
  // One entry = one run
  // ============================================================

  m_tree_runrate = new TTree(
      "RunRate",
      "Run-level MBD trigger and duration information");

  m_tree_runrate->Branch("run", &m_run, "run/I");

  m_tree_runrate->Branch(
      "run_bco_count", &m_run_bco_count, "run_bco_count/l");
  m_tree_runrate->Branch(
      "run_time_sec", &m_run_time_sec, "run_time_sec/D");

  m_tree_runrate->Branch("mbd_raw", &m_mbd_raw, "mbd_raw/l");
  m_tree_runrate->Branch("mbd_live", &m_mbd_live, "mbd_live/l");
  m_tree_runrate->Branch("mbd_scaled", &m_mbd_scaled, "mbd_scaled/l");

  m_tree_runrate->Branch(
      "mbd_scale_down", &m_mbd_scale_down, "mbd_scale_down/D");

  // ============================================================
  // Reset state
  // ============================================================

  m_first_event = true;
  m_bin_begin_bco = 0;
  m_last_bco = 0;
  m_segment = 0;

  ResetBin();

  m_run = -1;
  m_run_bco_count = 0;
  m_run_time_sec = 0.0;

  m_mbd_raw = 0;
  m_mbd_live = 0;
  m_mbd_scaled = 0;
  m_mbd_scale_down = 0.0;

  // ============================================================
  // TriggerRunInfo
  // ============================================================

  TriggerRunInfo *triggerRunInfo =
      findNode::getClass<TriggerRunInfo>(
          topNode, "TriggerRunInfo");

  if (!triggerRunInfo)
  {
    std::cout
        << "WARNING: TriggerRunInfo node not found"
        << std::endl;
  }
  else
  {
    // Clock scaler
    const int clock_bit =
        triggerRunInfo->getTriggerBitByName("Clock");

    if (clock_bit >= 0)
    {
      m_run_bco_count =
          static_cast<uint64_t>(
              triggerRunInfo->getRawScalersByBit(clock_bit));

      m_run_time_sec =
          static_cast<double>(m_run_bco_count) * BCO_SEC;
    }

    // MBD trigger
    const int mbd_bit =
        triggerRunInfo->getTriggerBitByName(
            "MBD N&S >= 1, vtx < 10 cm");

    if (mbd_bit >= 0)
    {
      m_mbd_raw =
          static_cast<uint64_t>(
              triggerRunInfo->getRawScalersByBit(mbd_bit));

      m_mbd_live =
          static_cast<uint64_t>(
              triggerRunInfo->getLiveScalersByBit(mbd_bit));

      m_mbd_scaled =
          static_cast<uint64_t>(
              triggerRunInfo->getScalersByBit(mbd_bit));

      m_mbd_scale_down =
          triggerRunInfo->getPrescaleByBit(mbd_bit);
    }
  }

  std::cout
      << "InttWordShiftCDBProduction settings:"
      << "\n  output = " << m_outfile_name
      << "\n  bin width [sec] = " << m_bin_width_sec
      << "\n  run time [sec] = " << m_run_time_sec
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttWordShiftCDBProduction::process_event(
    [[maybe_unused]] PHCompositeNode *topNode)
{
  const uint64_t bco =
      static_cast<uint64_t>(m_GL1Packet->getBCO());

  m_last_bco = bco;

  // First event
  if (m_first_event)
  {
    m_bin_begin_bco = bco;
    m_first_event = false;
  }

  // Check time bin
  const double elapsed_sec =
      static_cast<double>(bco - m_bin_begin_bco) * BCO_SEC;

  if (elapsed_sec >= m_bin_width_sec)
  {
    FinishBin(bco);
    m_bin_begin_bco = bco;
  }

  // ============================================================
  // Accumulate INTT clusters
  // ============================================================

  for (const auto &hitsetkey :
       m_TrkrClusterContainer->getHitSetKeys(
           TrkrDefs::TrkrId::inttId))
  {
    const int layer = TrkrDefs::getLayer(hitsetkey);

    if (layer < 3 || layer > 6)
    {
      continue;
    }

    const int phi = InttDefs::getLadderPhiId(hitsetkey);
    const int nphi = (layer < 5) ? 12 : 16;

    if (phi < 0 || phi >= nphi)
    {
      continue;
    }

    const int ladder_z = InttDefs::getLadderZId(hitsetkey);
    const int sn = (ladder_z < 2) ? 0 : 1;

    const int tb = InttDefs::getTimeBucketId(hitsetkey);

    auto clusrange =
        m_TrkrClusterContainer->getClusters(hitsetkey);

    for (auto iter = clusrange.first;
         iter != clusrange.second;
         ++iter)
    {
      TrkrCluster *cluster = iter->second;

      if (!cluster)
      {
        continue;
      }

      const double phi_size =
          static_cast<double>(cluster->getPhiSize());

      // TB = any
      ++m_all[layer][phi][sn].cluster_count;
      m_all[layer][phi][sn].phi_size_sum += phi_size;

      // TB = 0
      if (tb == 0)
      {
        ++m_tb0[layer][phi][sn].cluster_count;
        m_tb0[layer][phi][sn].phi_size_sum += phi_size;
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void InttWordShiftCDBProduction::FinishBin(
    const uint64_t end_bco)
{
  if (end_bco <= m_bin_begin_bco)
  {
    return;
  }

  const double bin_sec =
      static_cast<double>(end_bco - m_bin_begin_bco)
      * BCO_SEC;

  if (bin_sec <= 0.0)
  {
    return;
  }

  m_begin_bco = m_bin_begin_bco;
  m_end_bco = end_bco;

  // ============================================================
  // Save 112 half-ladders per bin
  // ============================================================

  for (int layer = 3; layer <= 6; ++layer)
  {
    const int nphi = (layer < 5) ? 12 : 16;

    for (int phi = 0; phi < nphi; ++phi)
    {
      for (int sn = 0; sn < MAX_SN; ++sn)
      {
        m_layer = layer;
        m_phi = phi;
        m_sn = sn;

        // ======================================================
        // FELIX mapping
        // ======================================================

        InttNameSpace::Offline_s off{};

        off.layer = layer;
        off.ladder_phi = phi;
        off.ladder_z = (sn == 0) ? 0 : 2;
        off.strip_x = 1;
        off.strip_y = 1;

        const InttNameSpace::RawData_s raw =
            InttNameSpace::ToRawData(off);

        m_fs = raw.felix_server;
        m_fc = raw.felix_channel;

        // ======================================================
        // TB = any
        // ======================================================

        const BinData &all = m_all[layer][phi][sn];

        m_cluster_count = all.cluster_count;
        m_cluster_rate = 0.0;
        m_mean_phi_size = 0.0;
        m_hit_rate = 0.0;

        if (all.cluster_count > 0)
        {
          m_cluster_rate =
              static_cast<double>(all.cluster_count)
              / bin_sec;

          m_mean_phi_size =
              all.phi_size_sum
              / static_cast<double>(all.cluster_count);

          m_hit_rate =
              m_cluster_rate * m_mean_phi_size;
        }

        // ======================================================
        // TB = 0
        // ======================================================

        const BinData &tb0 = m_tb0[layer][phi][sn];

        m_cluster_count_tb0 = tb0.cluster_count;
        m_cluster_rate_tb0 = 0.0;
        m_mean_phi_size_tb0 = 0.0;
        m_hit_rate_tb0 = 0.0;

        if (tb0.cluster_count > 0)
        {
          m_cluster_rate_tb0 =
              static_cast<double>(tb0.cluster_count)
              / bin_sec;

          m_mean_phi_size_tb0 =
              tb0.phi_size_sum
              / static_cast<double>(tb0.cluster_count);

          m_hit_rate_tb0 =
              m_cluster_rate_tb0 * m_mean_phi_size_tb0;
        }

        // Save statistics without bad/good decision
        m_tree_single->Fill();
      }
    }
  }

  ++m_segment;
  ResetBin();
}

//____________________________________________________________________________..
void InttWordShiftCDBProduction::ResetBin()
{
  for (int layer = 0; layer < MAX_LAYER; ++layer)
  {
    for (int phi = 0; phi < MAX_PHI; ++phi)
    {
      for (int sn = 0; sn < MAX_SN; ++sn)
      {
        m_all[layer][phi][sn].Reset();
        m_tb0[layer][phi][sn].Reset();
      }
    }
  }
}

//____________________________________________________________________________..
int InttWordShiftCDBProduction::ResetEvent(
    [[maybe_unused]] PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttWordShiftCDBProduction::EndRun(
    const int runnumber)
{
  m_run = runnumber;

  if (m_tree_runrate)
  {
    m_tree_runrate->Fill();
  }

  std::cout
      << "InttWordShiftCDBProduction::EndRun() run = "
      << runnumber
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttWordShiftCDBProduction::End(
    [[maybe_unused]] PHCompositeNode *topNode)
{
  // Save final partial bin
  if (!m_first_event && m_last_bco > m_bin_begin_bco)
  {
    FinishBin(m_last_bco);
  }

  if (!m_outfile)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_outfile->cd();

  if (m_tree_single)
  {
    m_tree_single->Write();
  }

  if (m_tree_runrate)
  {
    m_tree_runrate->Write();
  }

  m_outfile->Write();
  m_outfile->Close();

  std::cout
      << "Saved " << m_outfile_name
      << "\n  Single entries = "
      << (m_tree_single ? m_tree_single->GetEntries() : 0)
      << "\n  RunRate entries = "
      << (m_tree_runrate ? m_tree_runrate->GetEntries() : 0)
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttWordShiftCDBProduction::Reset(
    [[maybe_unused]] PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void InttWordShiftCDBProduction::Print(
    const std::string &what) const
{
  std::cout
      << "InttWordShiftCDBProduction::Print(): " << what
      << "\n  output = " << m_outfile_name
      << "\n  bin width [sec] = " << m_bin_width_sec
      << std::endl;
}

//____________________________________________________________________________..
void InttWordShiftCDBProduction::setOutFileName(
    const std::string &name)
{
  m_outfile_name = name;
}

//____________________________________________________________________________..
void InttWordShiftCDBProduction::setBinWidthSec(
    const double sec)
{
  m_bin_width_sec = sec;
}