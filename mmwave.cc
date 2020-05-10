/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/* *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Michele Polese <michele.polese@gmail.com>
 */

#include "ns3/mmwave-helper.h"
#include "ns3/epc-helper.h"
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"
#include "ns3/point-to-point-helper.h"
#include "ns3/config-store.h"




#include "ns3/mmwave-point-to-point-epc-helper.h"
//#include "ns3/gtk-config-store.h"
#include <ns3/buildings-helper.h>
#include <ns3/buildings-module.h>
#include <ns3/random-variable-stream.h>
#include <ns3/lte-ue-net-device.h>
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <list>


using namespace ns3;
using namespace mmwave;
using namespace std;

Ptr<PacketSink> sink;
uint64_t lastTotalRx = 0;
ofstream myFile;
NS_LOG_COMPONENT_DEFINE ("mmWaveTCPExample");
//Application to be installed on the UE
class MyAppTag : public Tag
{
public:
    MyAppTag ()
    {
    }

    MyAppTag (Time sendTs) : m_sendTs (sendTs)
    {
    }

    static TypeId GetTypeId (void)
    {
        static TypeId tid = TypeId ("ns3::MyAppTag")
                .SetParent<Tag> ()
                .AddConstructor<MyAppTag> ();
        return tid;
    }

    virtual TypeId  GetInstanceTypeId (void) const
    {
        return GetTypeId ();
    }

    virtual void  Serialize (TagBuffer i) const
    {
        i.WriteU64 (m_sendTs.GetNanoSeconds ());
    }

    virtual void  Deserialize (TagBuffer i)
    {
        m_sendTs = NanoSeconds (i.ReadU64 ());
    }

    virtual uint32_t  GetSerializedSize () const
    {
        return sizeof (m_sendTs);
    }

    virtual void Print (std::ostream &os) const
    {
        std::cout << m_sendTs;
    }

    Time m_sendTs;
};

class MyApp : public Application
{
public:
    MyApp ();
    virtual ~MyApp ();
    void ChangeDataRate (DataRate rate);
    void Setup (Ptr<Socket> socket, Address address, uint32_t packetSize, uint32_t nPackets, DataRate dataRate);



private:
    virtual void StartApplication (void);
    virtual void StopApplication (void);

    void ScheduleTx (void);
    void SendPacket (void);

    Ptr<Socket>     m_socket;
    Address         m_peer;
    uint32_t        m_packetSize;
    uint32_t        m_nPackets;
    DataRate        m_dataRate;
    EventId         m_sendEvent;
    bool            m_running;
    uint32_t        m_packetsSent;
};

MyApp::MyApp ()
        : m_socket (0),
          m_peer (),
          m_packetSize (0),
          m_nPackets (0),
          m_dataRate (0),
          m_sendEvent (),
          m_running (false),
          m_packetsSent (0)
{
}

MyApp::~MyApp ()
{
    m_socket = 0;
}

void
MyApp::Setup (Ptr<Socket> socket, Address address, uint32_t packetSize, uint32_t nPackets, DataRate dataRate)
{
    m_socket = socket;
    m_peer = address;
    m_packetSize = packetSize;
    m_nPackets = nPackets;
    m_dataRate = dataRate;
}

void
MyApp::ChangeDataRate (DataRate rate)
{
    m_dataRate = rate;
}

void
MyApp::StartApplication (void)
{
    m_running = true;
    m_packetsSent = 0;
    m_socket->Bind ();
    m_socket->Connect (m_peer);
    SendPacket ();
}

void
MyApp::StopApplication (void)
{
    m_running = false;

    if (m_sendEvent.IsRunning ())
    {
        Simulator::Cancel (m_sendEvent);
    }

    if (m_socket)
    {
        m_socket->Close ();
    }
}

void
MyApp::SendPacket (void)
{
    Ptr<Packet> packet = Create<Packet> (m_packetSize);
    MyAppTag tag (Simulator::Now ());

    m_socket->Send (packet);
    if (++m_packetsSent < m_nPackets)
    {
        ScheduleTx ();
    }
}



void
MyApp::ScheduleTx (void)
{
    if (m_running)
    {
        Time tNext (Seconds (m_packetSize * 8 / static_cast<double> (m_dataRate.GetBitRate ())));
        m_sendEvent = Simulator::Schedule (tNext, &MyApp::SendPacket, this);
    }
}


void
CalculateThroughput ()
{
    Time now = Simulator::Now ();                                         /* Return the simulator's virtual time. */
    double cur = (sink->GetTotalRx () - lastTotalRx) * (double) 8 / 1e5;     /* Convert Application RX Packets to MBits. */
    std::cout << now.GetSeconds () << "s: \t" << cur << " Mbit/s" << std::endl;
    myFile << now.GetSeconds () << "," << cur << std::endl;
    lastTotalRx = sink->GetTotalRx ();
    Simulator::Schedule (MilliSeconds (100), &CalculateThroughput);
}
//
//static void
//CwndChange (Ptr<OutputStreamWrapper> stream, uint32_t oldCwnd, uint32_t newCwnd)
//{
//    *stream->GetStream () << Simulator::Now ().GetSeconds () << "\t" << oldCwnd << "\t" << newCwnd << std::endl;
//}
//
//static void
//RttChange (Ptr<OutputStreamWrapper> stream, Time oldRtt, Time newRtt)
//{
//    *stream->GetStream () << Simulator::Now ().GetSeconds () << "\t" << oldRtt.GetSeconds () << "\t" << newRtt.GetSeconds () << std::endl;
//}


//static void
//RxDrop (Ptr<PcapFileWrapper> file, Ptr<const Packet> p)
//{
//    file->Write (Simulator::Now (), p);
//    NS_LOG_UNCOND (p);
//
//}
//
//static void Rx (Ptr<OutputStreamWrapper> stream, Ptr<const Packet> packet, const Address &from)
//{
//    *stream->GetStream () << Simulator::Now ().GetSeconds () << "\t" << packet->GetSize () << std::endl;
//}

void CourseChange(std::string context, Ptr<const MobilityModel> model)
{
    Vector position = model->GetPosition ();
    NS_LOG_UNCOND (context << " x = " << position.x << ",y = " << position.y);
}

void
PrintPosition (Ptr<Node> node)
{
    Ptr<MobilityModel> model = node->GetObject<MobilityModel> ();
    NS_LOG_UNCOND ("Position +****************************** " << model->GetPosition () << " at time " << Simulator::Now ().GetSeconds ());
}

void
ChangeSpeed (Ptr<Node> n, Vector speed)
{
    n->GetObject<ConstantVelocityMobilityModel> ()->SetVelocity (speed);
    NS_LOG_UNCOND ("************************--------------------Change Speed-------------------------------*****************");
}


void
PrintDistance (Ptr<Node> node, NodeContainer enbNodes)
{
    Ptr<MobilityModel> model = node->GetObject<MobilityModel> ();
    for (uint32_t i = 0; i < enbNodes.GetN() ; i++)
        NS_LOG_UNCOND ("Distance" << i << " +****************************** " << CalculateDistance(model->GetPosition (), enbNodes.Get (i)->GetObject<MobilityModel>()->GetPosition()) << " at time " << Simulator::Now ().GetSeconds ());
}

bool
AreOverlapping (Box a, Box b)
{
    return !((a.xMin > b.xMax) || (b.xMin > a.xMax) || (a.yMin > b.yMax) || (b.yMin > a.yMax) );
}


bool
OverlapWithAnyPrevious (Box box, std::list<Box> m_previousBlocks)
{
    for (std::list<Box>::iterator it = m_previousBlocks.begin (); it != m_previousBlocks.end (); ++it)
    {
        if (AreOverlapping (*it,box))
        {
            return true;
        }
    }
    return false;
}

std::pair<Box, std::list<Box> >
GenerateBuildingBounds (double xArea, double yArea, double maxBuildSize, std::list<Box> m_previousBlocks )
{

    Ptr<UniformRandomVariable> xMinBuilding = CreateObject<UniformRandomVariable> ();
    xMinBuilding->SetAttribute ("Min",DoubleValue (30));
    xMinBuilding->SetAttribute ("Max",DoubleValue (xArea));

    NS_LOG_UNCOND ("min " << 0 << " max " << xArea);

    Ptr<UniformRandomVariable> yMinBuilding = CreateObject<UniformRandomVariable> ();
    yMinBuilding->SetAttribute ("Min",DoubleValue (0));
    yMinBuilding->SetAttribute ("Max",DoubleValue (yArea));

    NS_LOG_UNCOND ("min " << 0 << " max " << yArea);

    Box box;
    uint32_t attempt = 0;
    do
    {
        NS_ASSERT_MSG (attempt < 100, "Too many failed attempts to position non-overlapping buildings. Maybe area too small or too many buildings?");
        box.xMin = xMinBuilding->GetValue ();

        Ptr<UniformRandomVariable> xMaxBuilding = CreateObject<UniformRandomVariable> ();
        xMaxBuilding->SetAttribute ("Min",DoubleValue (box.xMin));
        xMaxBuilding->SetAttribute ("Max",DoubleValue (box.xMin + maxBuildSize));
        box.xMax = xMaxBuilding->GetValue ();

        box.yMin = yMinBuilding->GetValue ();

        Ptr<UniformRandomVariable> yMaxBuilding = CreateObject<UniformRandomVariable> ();
        yMaxBuilding->SetAttribute ("Min",DoubleValue (box.yMin));
        yMaxBuilding->SetAttribute ("Max",DoubleValue (box.yMin + maxBuildSize));
        box.yMax = yMaxBuilding->GetValue ();

        ++attempt;
    }
    while (OverlapWithAnyPrevious (box, m_previousBlocks));


    NS_LOG_UNCOND ("Building in coordinates (" << box.xMin << " , " << box.yMin << ") and ("  << box.xMax << " , " << box.yMax <<
                                               ") accepted after " << attempt << " attempts");
    m_previousBlocks.push_back (box);
    std::pair<Box, std::list<Box> > pairReturn = std::make_pair (box,m_previousBlocks);
    return pairReturn;

}


/**
PLOTTING
**/
void
PrintGnuplottableBuildingListToFile (std::string filename)
{
    std::ofstream outFile;
    outFile.open (filename.c_str (), std::ios_base::out | std::ios_base::trunc);
    if (!outFile.is_open ())
    {
        NS_LOG_ERROR ("Can't open file " << filename);
        return;
    }
    uint32_t index = 0;
    for (BuildingList::Iterator it = BuildingList::Begin (); it != BuildingList::End (); ++it)
    {
        ++index;
        Box box = (*it)->GetBoundaries ();
        outFile << "set object " << index
                << " rect from " << box.xMin  << "," << box.yMin
                << " to "   << box.xMax  << "," << box.yMax
                << " front fs empty "
                << std::endl;
    }
}

void
PrintGnuplottableUeListToFile (std::string filename)
{
    std::ofstream outFile;
    outFile.open (filename.c_str (), std::ios_base::out | std::ios_base::trunc);
    if (!outFile.is_open ())
    {
        NS_LOG_ERROR ("Can't open file " << filename);
        return;
    }
    for (NodeList::Iterator it = NodeList::Begin (); it != NodeList::End (); ++it)
    {
        Ptr<Node> node = *it;
        int nDevs = node->GetNDevices ();
        for (int j = 0; j < nDevs; j++)
        {
            Ptr<LteUeNetDevice> uedev = node->GetDevice (j)->GetObject <LteUeNetDevice> ();
            Ptr<MmWaveUeNetDevice> mmuedev = node->GetDevice (j)->GetObject <MmWaveUeNetDevice> ();
            Ptr<McUeNetDevice> mcuedev = node->GetDevice (j)->GetObject <McUeNetDevice> ();
            // NS_LOG_UNCOND (node->GetDevice (j));

            if (uedev)
            {
                Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
                outFile << "set label \"" << uedev->GetImsi ()
                        << "\" at " << pos.x << "," << pos.y << " left font \"Helvetica,8\" textcolor rgb \"black\" front point pt 1 ps 0.3 lc rgb \"black\" offset 0,0"
                        << std::endl;
            }
            else if (mmuedev)
            {
                Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
                outFile << "set label \"" << mmuedev->GetImsi ()
                        << "\" at " << pos.x << "," << pos.y << " left font \"Helvetica,8\" textcolor rgb \"black\" front point pt 1 ps 0.3 lc rgb \"black\" offset 0,0"
                        << std::endl;
            }
            else if (mcuedev)
            {
                Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
                outFile << "set label \"" << mcuedev->GetImsi ()
                        << "\" at " << pos.x << "," << pos.y << " left font \"Helvetica,8\" textcolor rgb \"black\" front point pt 1 ps 0.3 lc rgb \"black\" offset 0,0"
                        << std::endl;
            }
        }
    }
}

void
PrintGnuplottableEnbListToFile (std::string filename)
{
    std::ofstream outFile;
    outFile.open (filename.c_str (), std::ios_base::out | std::ios_base::trunc);
    if (!outFile.is_open ())
    {
        NS_LOG_ERROR ("Can't open file " << filename);
        return;
    }
    for (NodeList::Iterator it = NodeList::Begin (); it != NodeList::End (); ++it)
    {
        Ptr<Node> node = *it;

        int nDevs = node->GetNDevices ();
        NS_LOG_UNCOND(nDevs);
        for (int j = 0; j < nDevs; j++)
        {
            // Ptr<MmWaveEnbNetDevice> mm2 = CreateObject<MmWaveEnbNetDevice> ();
            Ptr<LteEnbNetDevice> enbdev = node->GetDevice (j)->GetObject <LteEnbNetDevice> ();
            Ptr<MmWaveEnbNetDevice> mmdev = node->GetDevice (j)->GetObject <MmWaveEnbNetDevice> ();
            if (enbdev)
            {
                Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
                outFile << "set label \"" << enbdev->GetCellId ()
                        << "\" at " << pos.x << "," << pos.y
                        << " left font \"Helvetica,8\" textcolor rgb \"blue\" front  point pt 4 ps 0.3 lc rgb \"blue\" offset 0,0"
                        << std::endl;
            }
            else if (mmdev)
            {
                Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
                outFile << "set label \"" << mmdev->GetCellId ()
                        << "\" at " << pos.x << "," << pos.y
                        << " left font \"Helvetica,8\" textcolor rgb \"red\" front  point pt 4 ps 0.3 lc rgb \"red\" offset 0,0"
                        << std::endl;
            }

        }
    }
}


static ns3::GlobalValue g_mmw1DistFromMainStreet ("mmw1Dist", "Distance from the main street of the first MmWaveEnb",
                                                  ns3::UintegerValue (50), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_mmw2DistFromMainStreet ("mmw2Dist", "Distance from the main street of the second MmWaveEnb",
                                                  ns3::UintegerValue (50), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_mmw3DistFromMainStreet ("mmw3Dist", "Distance from the main street of the third MmWaveEnb",
                                                  ns3::UintegerValue (110), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_mmWaveDistance ("mmWaveDist", "Distance between MmWave eNB 1 and 2",
                                          ns3::UintegerValue (200), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_interPckInterval ("interPckInterval", "Interarrival time of UDP packets (us)",
                                            ns3::UintegerValue (20), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_bufferSize ("bufferSize", "RLC tx buffer size (MB)",
                                      ns3::UintegerValue (20), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_x2Latency ("x2Latency", "Latency on X2 interface (us)",
                                     ns3::DoubleValue (500), ns3::MakeDoubleChecker<double> ());
static ns3::GlobalValue g_mmeLatency ("mmeLatency", "Latency on MME interface (us)",
                                      ns3::DoubleValue (10000), ns3::MakeDoubleChecker<double> ());
static ns3::GlobalValue g_mobileUeSpeed ("mobileSpeed", "The speed of the UE (m/s)",
                                         ns3::DoubleValue (20), ns3::MakeDoubleChecker<double> ());
static ns3::GlobalValue g_rlcAmEnabled ("rlcAmEnabled", "If true, use RLC AM, else use RLC UM",
                                        ns3::BooleanValue (true), ns3::MakeBooleanChecker ());
static ns3::GlobalValue g_runNumber ("runNumber", "Run number for rng",
                                     ns3::UintegerValue (10), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_outPath ("outPath",
                                   "The path of output log files",
                                   ns3::StringValue ("./"), ns3::MakeStringChecker ());
static ns3::GlobalValue g_noiseAndFilter ("noiseAndFilter", "If true, use noisy SINR samples, filtered. If false, just use the SINR measure",
                                          ns3::BooleanValue (false), ns3::MakeBooleanChecker ());
static ns3::GlobalValue g_handoverMode ("handoverMode",
                                        "Handover mode",
                                        ns3::UintegerValue (2), ns3::MakeUintegerChecker<uint8_t> ());
static ns3::GlobalValue g_reportTablePeriodicity ("reportTablePeriodicity", "Periodicity of RTs",
                                                  ns3::UintegerValue (1600), ns3::MakeUintegerChecker<uint32_t> ());
static ns3::GlobalValue g_outageThreshold ("outageTh", "Outage threshold",
                                           ns3::DoubleValue (-5), ns3::MakeDoubleChecker<double> ());
static ns3::GlobalValue g_lteUplink ("lteUplink", "If true, always use LTE for uplink signalling",
                                     ns3::BooleanValue (false), ns3::MakeBooleanChecker ());

int
main (int argc, char *argv[])
{

  bool harqEnabled = true;
  bool fixedTti = false;
  unsigned symPerSf = 24;
  double sfPeriod = 100.0;

  std::list<Box>  m_previousBlocks;

  // Command line arguments
  CommandLine cmd;
  cmd.Parse (argc, argv);

  UintegerValue uintegerValue;
  BooleanValue booleanValue;
  StringValue stringValue;
  DoubleValue doubleValue;

  double ueInitialPosition = -200;
  double ueFinalPosition = 600;

  // Variables for the RT
  int windowForTransient = 150; // number of samples for the vector to use in the filter
  GlobalValue::GetValueByName ("reportTablePeriodicity", uintegerValue);
  int ReportTablePeriodicity = (int)uintegerValue.Get (); // in microseconds
  if (ReportTablePeriodicity == 1600)
    {
      windowForTransient = 150;
    }
  else if (ReportTablePeriodicity == 25600)
    {
      windowForTransient = 50;
    }
  else if (ReportTablePeriodicity == 12800)
    {
      windowForTransient = 100;
    }
  else
    {
      NS_ASSERT_MSG (false, "Unrecognized");
    }

  int vectorTransient = windowForTransient * ReportTablePeriodicity;

  // params for RT, filter, HO mode
  GlobalValue::GetValueByName ("noiseAndFilter", booleanValue);
  bool noiseAndFilter = booleanValue.Get ();
  GlobalValue::GetValueByName ("handoverMode", uintegerValue);
  uint8_t hoMode = uintegerValue.Get ();
  GlobalValue::GetValueByName ("outageTh", doubleValue);
  double outageTh = doubleValue.Get ();

  GlobalValue::GetValueByName ("rlcAmEnabled", booleanValue);
  bool rlcAmEnabled = booleanValue.Get ();
  GlobalValue::GetValueByName ("bufferSize", uintegerValue);
  uint32_t bufferSize = uintegerValue.Get ();
  GlobalValue::GetValueByName ("interPckInterval", uintegerValue);
  uint32_t interPacketInterval = uintegerValue.Get ();
  GlobalValue::GetValueByName ("x2Latency", doubleValue);
  double x2Latency = doubleValue.Get ();
  GlobalValue::GetValueByName ("mmeLatency", doubleValue);
  double mmeLatency = doubleValue.Get ();
  GlobalValue::GetValueByName ("mobileSpeed", doubleValue);
  double ueSpeed = doubleValue.Get ();

  double transientDuration = double(vectorTransient) / 1000000;
  double simTime = transientDuration + ((double)ueFinalPosition - (double)ueInitialPosition) / ueSpeed + 1;

  // rng things
  GlobalValue::GetValueByName ("runNumber", uintegerValue);
  uint32_t runSet = uintegerValue.Get ();
  uint32_t seedSet = 5;
  RngSeedManager::SetSeed (seedSet);
  RngSeedManager::SetRun (runSet);
  char seedSetStr[21];
  char runSetStr[21];
  sprintf (seedSetStr, "%d", seedSet);
  sprintf (runSetStr, "%d", runSet);

  GlobalValue::GetValueByName ("outPath", stringValue);
  std::string path = stringValue.Get ();
//  std::string mmWaveOutName = "MmWaveSwitchStatsX";
//  std::string lteOutName = "LteSwitchStatsX";
//  std::string dlRlcOutName = "DlRlcStatsX";
//  std::string dlPdcpOutName = "DlPdcpStatsX";
//  std::string ulRlcOutName = "UlRlcStatsX";
//  std::string ulPdcpOutName = "UlPdcpStatsX";
  std::string  ueHandoverStartOutName =  "UeHandoverStartStats";
//  std::string enbHandoverStartOutName = "EnbHandoverStartStats";
  std::string  ueHandoverEndOutName =  "UeHandoverEndStats";
//  std::string enbHandoverEndOutName = "EnbHandoverEndStatsX";
//  std::string cellIdInTimeOutName = "CellIdStatsX";
//  std::string cellIdInTimeHandoverOutName = "CellIdStatsHandoverX";
 std::string mmWaveSinrOutputFilename = "MmWaveSinrTime";
//  std::string x2statOutputFilename = "X2StatsX";
//  std::string udpSentFilename = "UdpSentX";
//  std::string udpReceivedFilename = "UdpReceivedX";
  std::string extension = ".txt";
    std::string version;
    std::string handoff;
    switch (hoMode)
    {
        case 1:
            handoff = "threshold";
            break;
        case 2:
            handoff = "fix";;
            break;
        case 3:
            handoff = "dyn";
            break;
    }

    version = "blarge3bsv3ttt50" + handoff + std::to_string((int)ueSpeed);
    myFile.open(version + "MmWaveThroughput.txt");


    NS_LOG_UNCOND ("rlcAmEnabled " << rlcAmEnabled << " bufferSize " << bufferSize << " interPacketInterval " <<
                                   interPacketInterval << " x2Latency " << x2Latency << " mmeLatency " << mmeLatency << " mobileSpeed " << ueSpeed << " fileName " << version);

  Config::SetDefault ("ns3::MmWaveUeMac::UpdateUeSinrEstimatePeriod", DoubleValue (0));

  //get current time
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  strftime (buffer,80,"%d_%m_%Y_%I_%M_%S",timeinfo);
  std::string time_str (buffer);

  Config::SetDefault ("ns3::MmWaveHelper::RlcAmEnabled", BooleanValue (rlcAmEnabled));
  Config::SetDefault ("ns3::MmWaveHelper::HarqEnabled", BooleanValue (harqEnabled));
  Config::SetDefault ("ns3::MmWaveFlexTtiMacScheduler::HarqEnabled", BooleanValue (harqEnabled));
  Config::SetDefault ("ns3::MmWaveFlexTtiMaxWeightMacScheduler::HarqEnabled", BooleanValue (harqEnabled));
  Config::SetDefault ("ns3::MmWaveFlexTtiMaxWeightMacScheduler::FixedTti", BooleanValue (fixedTti));
  Config::SetDefault ("ns3::MmWaveFlexTtiMaxWeightMacScheduler::SymPerSlot", UintegerValue (6));
  Config::SetDefault ("ns3::MmWavePhyMacCommon::ResourceBlockNum", UintegerValue (1));
  Config::SetDefault ("ns3::MmWavePhyMacCommon::ChunkPerRB", UintegerValue (72));
  Config::SetDefault ("ns3::MmWavePhyMacCommon::SymbolsPerSubframe", UintegerValue (symPerSf));
  Config::SetDefault ("ns3::MmWavePhyMacCommon::SubframePeriod", DoubleValue (sfPeriod));
  Config::SetDefault ("ns3::MmWavePhyMacCommon::TbDecodeLatency", UintegerValue (200.0));
  Config::SetDefault ("ns3::MmWavePhyMacCommon::NumHarqProcess", UintegerValue (100));
  Config::SetDefault ("ns3::MmWaveBeamforming::LongTermUpdatePeriod", TimeValue (MilliSeconds (100.0)));
  Config::SetDefault ("ns3::LteEnbRrc::SystemInformationPeriodicity", TimeValue (MilliSeconds (5.0)));
  Config::SetDefault ("ns3::LteRlcAm::ReportBufferStatusTimer", TimeValue (MicroSeconds (100.0)));
  Config::SetDefault ("ns3::LteRlcUmLowLat::ReportBufferStatusTimer", TimeValue (MicroSeconds (100.0)));
  Config::SetDefault ("ns3::LteEnbRrc::SrsPeriodicity", UintegerValue (320));
  Config::SetDefault ("ns3::LteEnbRrc::FirstSibTime", UintegerValue (2));
  Config::SetDefault ("ns3::MmWavePointToPointEpcHelper::X2LinkDelay", TimeValue (MicroSeconds (x2Latency)));
  Config::SetDefault ("ns3::MmWavePointToPointEpcHelper::X2LinkDataRate", DataRateValue (DataRate ("1000Gb/s")));
  Config::SetDefault ("ns3::MmWavePointToPointEpcHelper::X2LinkMtu",  UintegerValue (10000));
  Config::SetDefault ("ns3::MmWavePointToPointEpcHelper::S1uLinkDelay", TimeValue (MicroSeconds (1000)));
  Config::SetDefault ("ns3::MmWavePointToPointEpcHelper::S1apLinkDelay", TimeValue (MicroSeconds (mmeLatency)));
//  Config::SetDefault ("ns3::McStatsCalculator::MmWaveOutputFilename", StringValue                 (path + version + mmWaveOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::McStatsCalculator::LteOutputFilename", StringValue                    (path + version + lteOutName    + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::McStatsCalculator::CellIdInTimeOutputFilename", StringValue           (path + version + cellIdInTimeOutName    + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::MmWaveBearerStatsCalculator::DlRlcOutputFilename", StringValue        (path + version + dlRlcOutName   + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  //Config::SetDefault ("ns3::MmWaveBearerStmatsCalculator::UlRlcOutputFilename", StringValue        (path + version + ulRlcOutName   + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::MmWaveBearerStatsCalculator::DlPdcpOutputFilename", StringValue       (path + version + dlPdcpOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::MmWaveBearerStatsCalculator::UlPdcpOutputFilename", StringValue       (path + version + ulPdcpOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
Config::SetDefault ("ns3::MmWaveBearerStatsConnector::UeHandoverStartOutputFilename", StringValue    (path + version +  ueHandoverStartOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::MmWaveBearerStatsConnector::EnbHandoverStartOutputFilename", StringValue   (path + version + enbHandoverStartOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
Config::SetDefault ("ns3::MmWaveBearerStatsConnector::UeHandoverEndOutputFilename", StringValue    (path + version +  ueHandoverEndOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::MmWaveBearerStatsConnector::EnbHandoverEndOutputFilename", StringValue   (path + version + enbHandoverEndOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::MmWaveBearerStatsConnector::CellIdStatsHandoverOutputFilename", StringValue (path + version + cellIdInTimeHandoverOutName + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
Config::SetDefault ("ns3::MmWaveBearerStatsConnector::MmWaveSinrOutputFilename", StringValue (path + version + mmWaveSinrOutputFilename + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::CoreNetworkStatsCalculator::X2FileName", StringValue                  (path + version + x2statOutputFilename    + "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));

  //  std::string lostFilename = path + version + "LostUdpPackets" +  "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension;
//  Config::SetDefault ("ns3::UdpServer::ReceivedPacketsFilename", StringValue(path + version + "ReceivedUdp" +  "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::UdpClient::SentPacketsFilename", StringValue(path + version + "SentUdp" +  "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::UdpServer::ReceivedSnFilename", StringValue(path + version + "ReceivedSn" +  "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));
//  Config::SetDefault ("ns3::LteRlcAm::BufferSizeFilename", StringValue (path + version + "RlcAmBufferSize" +  "_" + seedSetStr + "_" + runSetStr + "_" + time_str + extension));

  Config::SetDefault ("ns3::LteRlcUm::MaxTxBufferSize", UintegerValue (bufferSize * 1024 * 1024));
  Config::SetDefault ("ns3::LteRlcUmLowLat::MaxTxBufferSize", UintegerValue (bufferSize * 1024 * 1024));
  Config::SetDefault ("ns3::LteRlcAm::StatusProhibitTimer", TimeValue (MilliSeconds (10.0)));
  Config::SetDefault ("ns3::LteRlcAm::MaxTxBufferSize", UintegerValue (bufferSize * 1024 * 1024));

  // handover and RT related params
  switch (hoMode)
    {
    case 1:
      Config::SetDefault ("ns3::LteEnbRrc::SecondaryCellHandoverMode", EnumValue (LteEnbRrc::THRESHOLD));
      break;
    case 2:
      Config::SetDefault ("ns3::LteEnbRrc::SecondaryCellHandoverMode", EnumValue (LteEnbRrc::FIXED_TTT));
      break;
    case 3:
      Config::SetDefault ("ns3::LteEnbRrc::SecondaryCellHandoverMode", EnumValue (LteEnbRrc::DYNAMIC_TTT));
      break;
    }

  Config::SetDefault ("ns3::LteEnbRrc::FixedTttValue", UintegerValue (50));
  Config::SetDefault ("ns3::LteEnbRrc::CrtPeriod", IntegerValue (ReportTablePeriodicity));
  Config::SetDefault ("ns3::LteEnbRrc::OutageThreshold", DoubleValue (outageTh));
  Config::SetDefault ("ns3::MmWaveEnbPhy::UpdateSinrEstimatePeriod", IntegerValue (ReportTablePeriodicity));
  Config::SetDefault ("ns3::MmWaveEnbPhy::Transient", IntegerValue (vectorTransient));
  Config::SetDefault ("ns3::MmWaveEnbPhy::NoiseAndFilter", BooleanValue (noiseAndFilter));

  GlobalValue::GetValueByName ("lteUplink", booleanValue);
  bool lteUplink = booleanValue.Get ();

  Config::SetDefault ("ns3::McUePdcp::LteUplink", BooleanValue (lteUplink));
  std::cout << "Lte uplink " << lteUplink << "\n";

  // settings for the 3GPP the channel
  Config::SetDefault ("ns3::MmWave3gppPropagationLossModel::ChannelCondition", StringValue ("a"));
  Config::SetDefault ("ns3::MmWave3gppPropagationLossModel::Scenario", StringValue ("UMa"));
  Config::SetDefault ("ns3::MmWave3gppPropagationLossModel::OptionalNlos", BooleanValue (true));
  Config::SetDefault ("ns3::MmWave3gppPropagationLossModel::Shadowing", BooleanValue (true)); // enable or disable the shadowing effect
  Config::SetDefault ("ns3::MmWave3gppBuildingsPropagationLossModel::UpdateCondition", BooleanValue (true)); // enable or disable the LOS/NLOS update when the UE move
  Config::SetDefault ("ns3::AntennaArrayModel::AntennaHorizontalSpacing", DoubleValue (0.5));
  Config::SetDefault ("ns3::AntennaArrayModel::AntennaVerticalSpacing", DoubleValue (0.5));
  Config::SetDefault ("ns3::MmWave3gppChannel::UpdatePeriod", TimeValue (MilliSeconds (100))); // interval after which the channel for a moving user is updated,
  // with spatial consistency procedure. If 0, spatial consistency is not used
  Config::SetDefault ("ns3::MmWave3gppChannel::DirectBeam", BooleanValue (true)); // Set true to perform the beam in the exact direction of receiver node.
  Config::SetDefault ("ns3::MmWave3gppChannel::Blockage", BooleanValue (true)); // use blockage or not
  Config::SetDefault ("ns3::MmWave3gppChannel::PortraitMode", BooleanValue (true)); // use blockage model with UT in portrait mode
  Config::SetDefault ("ns3::MmWave3gppChannel::NumNonselfBlocking", IntegerValue (4)); // number of non-self blocking obstacles

  // set the number of antennas in the devices
  Config::SetDefault ("ns3::McUeNetDevice::AntennaNum", UintegerValue(16));
  Config::SetDefault ("ns3::MmWaveEnbNetDevice::AntennaNum", UintegerValue(64));

  Ptr<MmWaveHelper> mmwaveHelper = CreateObject<MmWaveHelper> ();
  if (true)
    {
      mmwaveHelper->SetAttribute ("PathlossModel", StringValue ("ns3::MmWave3gppBuildingsPropagationLossModel"));

    }
  else
    {
      mmwaveHelper->SetAttribute ("PathlossModel", StringValue ("ns3::MmWave3gppPropagationLossModel"));
    }

  mmwaveHelper->SetAttribute ("ChannelModel", StringValue ("ns3::MmWave3gppChannel"));


  //Self Added
  Config::SetDefault ("ns3::MmWaveEnbPhy::TxPower", DoubleValue (30));
  Config::SetDefault ("ns3::MmWavePhyMacCommon::CenterFreq", DoubleValue (100e9));
  //Ptr<MmWaveHelper> mmwaveHelper = CreateObject<MmWaveHelper> ();
  //mmwaveHelper->SetSchedulerType ("ns3::MmWaveFlexTtiMaxWeightMacScheduler");
  Ptr<MmWavePointToPointEpcHelper> epcHelper = CreateObject<MmWavePointToPointEpcHelper> ();
  mmwaveHelper->SetEpcHelper (epcHelper);
  mmwaveHelper->SetHarqEnabled (harqEnabled);
//  mmwaveHelper->SetAttribute ("PathlossModel", StringValue ("ns3::BuildingsObstaclePropagationLossModel"));
  mmwaveHelper->Initialize ();

  ConfigStore inputConfig;
  inputConfig.ConfigureDefaults ();

  // parse again so you can override default values from the command line
  cmd.Parse (argc, argv);

  // Get SGW/PGW and create a single RemoteHost
  Ptr<Node> pgw = epcHelper->GetPgwNode ();
  NodeContainer remoteHostContainer;
  remoteHostContainer.Create (1);
  Ptr<Node> remoteHost = remoteHostContainer.Get (0);
  InternetStackHelper internet;
  internet.Install (remoteHostContainer);

  // Create the Internet by connecting remoteHost to pgw. Setup routing too
  PointToPointHelper p2ph;
  p2ph.SetDeviceAttribute ("DataRate", DataRateValue (DataRate ("100Gb/s")));
  p2ph.SetDeviceAttribute ("Mtu", UintegerValue (2500));
  p2ph.SetChannelAttribute ("Delay", TimeValue (Seconds (0.010)));
  NetDeviceContainer internetDevices = p2ph.Install (pgw, remoteHost);
  Ipv4AddressHelper ipv4h;
  ipv4h.SetBase ("1.0.0.0", "255.0.0.0");
  Ipv4InterfaceContainer internetIpIfaces = ipv4h.Assign (internetDevices);
  // interface 0 is localhost, 1 is the p2p device
  //Ipv4Address remoteHostAddr = internetIpIfaces.GetAddress (1);
  Ipv4StaticRoutingHelper ipv4RoutingHelper;
  Ptr<Ipv4StaticRouting> remoteHostStaticRouting = ipv4RoutingHelper.GetStaticRouting (remoteHost->GetObject<Ipv4> ());
  remoteHostStaticRouting->AddNetworkRouteTo (Ipv4Address ("7.0.0.0"), Ipv4Mask ("255.0.0.0"), 1);

  // create LTE, mmWave eNB nodes and UE node
  NodeContainer ueNodes;
  NodeContainer mmWaveEnbNodes;
  NodeContainer lteEnbNodes;
  NodeContainer allEnbNodes;
  mmWaveEnbNodes.Create (4);
  lteEnbNodes.Create (1);
  ueNodes.Create (1);
  allEnbNodes.Add (lteEnbNodes);
  allEnbNodes.Add (mmWaveEnbNodes);

  // Positions of the eNBs
  Vector ltePosition = Vector (20000, 70, 3);
  Vector mmw1Position = Vector (0, 70, 3); //was 50
 Vector mmw2Position = Vector (200, 70, 3); //was 150
    Vector mmw3Position = Vector (400, 70, 3);
Vector mmw4Position = Vector (-200, 70, 3); //was 50
//  Vector mmw5Position = Vector (80, 70, 3); //was 150
//  Vector mmw6Position = Vector (100, 70, 3);


//Position the Buildings
    Ptr < Building > building2;
    building2 = Create<Building> ();
    building2->SetBoundaries (Box (190.0, 210.0,
                                   100.0, 140.0,
                                   0.0, 20.0));
    building2->SetBuildingType (Building::Residential);
    building2->SetExtWallsType (Building::ConcreteWithoutWindows);
    building2->SetNFloors (1);
    building2->SetNRoomsX (1);
    building2->SetNRoomsY (1);

//    Ptr < Building > building3;
//    building3 = Create<Building> ();
//    building3->SetBoundaries (Box (380.0, 420.0,
//                                   40.0, 30.0,
//                                   0.0, 20.0));
//    building3->SetBuildingType (Building::Residential);
//    building3->SetExtWallsType (Building::ConcreteWithoutWindows);
//    building3->SetNFloors (1);
//    building3->SetNRoomsX (1);
//    building3->SetNRoomsY (1);

//    Ptr < Building > building4;
//    building4 = Create<Building> ();
//    building4->SetBoundaries (Box (580.0, 620.0,
//                                    115.0, 125.0,
//                                    0.0, 20.0));
//    building4->SetBuildingType (Building::Residential);
//    building4->SetExtWallsType (Building::ConcreteWithoutWindows);
//    building4->SetNFloors (1);
//    building4->SetNRoomsX (1);
//    building4->SetNRoomsY (1);
//
//
//
//
//    Ptr < Building > building5;
//    building5 = Create<Building> ();
//    building5->SetBoundaries (Box (88.0, 92.0,
//                                   68.0, 72.0,
//                                   0.0, 20.0));
//    building5->SetBuildingType (Building::Residential);
//    building5->SetExtWallsType (Building::ConcreteWithoutWindows);
//    building5->SetNFloors (1);
//    building5->SetNRoomsX (1);
//    building5->SetNRoomsY (1);
    // Install Mobility Model


    // Install Mobility Model
  Ptr<ListPositionAllocator> enbPositionAlloc = CreateObject<ListPositionAllocator> ();
  enbPositionAlloc->Add (ltePosition); // LTE BS, out of area where buildings are deployed
  enbPositionAlloc->Add (mmw1Position);
  enbPositionAlloc->Add (mmw2Position);
enbPositionAlloc->Add (mmw3Position);
enbPositionAlloc->Add (mmw4Position);
//  enbPositionAlloc->Add (mmw5Position);
//  enbPositionAlloc->Add (mmw6Position);
  MobilityHelper enbmobility;
  enbmobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  enbmobility.SetPositionAllocator (enbPositionAlloc);
  enbmobility.Install (allEnbNodes);
  BuildingsHelper::Install (allEnbNodes);

  MobilityHelper uemobility;
  Ptr<ListPositionAllocator> uePositionAlloc = CreateObject<ListPositionAllocator> ();
  //uePositionAlloc->Add (Vector (ueInitialPosition, -5, 0));
  uePositionAlloc->Add (Vector (ueInitialPosition, 170, 1.6));
  uemobility.SetMobilityModel ("ns3::ConstantVelocityMobilityModel");
  uemobility.SetPositionAllocator (uePositionAlloc);
  uemobility.Install (ueNodes);
  BuildingsHelper::Install (ueNodes);

  //ueNodes.Get (0)->GetObject<MobilityModel> ()->SetPosition (Vector (ueInitialPosition, -5, 0));
  ueNodes.Get (0)->GetObject<MobilityModel> ()->SetPosition (Vector (ueInitialPosition,170, 1.6));
  ueNodes.Get (0)->GetObject<ConstantVelocityMobilityModel> ()->SetVelocity (Vector (0, 0, 0));

  // Install mmWave, lte, mc Devices to the nodes
  NetDeviceContainer lteEnbDevs = mmwaveHelper->InstallLteEnbDevice (lteEnbNodes);
  NetDeviceContainer mmWaveEnbDevs = mmwaveHelper->InstallEnbDevice (mmWaveEnbNodes);
  NetDeviceContainer mcUeDevs;
  mcUeDevs = mmwaveHelper->InstallMcUeDevice (ueNodes);

  // Install the IP stack on the UEs
  internet.Install (ueNodes);
  Ipv4InterfaceContainer ueIpIface;
  ueIpIface = epcHelper->AssignUeIpv4Address (NetDeviceContainer (mcUeDevs));
  // Assign IP address to UEs, and install applications
  for (uint32_t u = 0; u < ueNodes.GetN (); ++u)
    {
      Ptr<Node> ueNode = ueNodes.Get (u);
      // Set the default gateway for the UE
      Ptr<Ipv4StaticRouting> ueStaticRouting = ipv4RoutingHelper.GetStaticRouting (ueNode->GetObject<Ipv4> ());
      ueStaticRouting->SetDefaultRoute (epcHelper->GetUeDefaultGatewayAddress (), 1);
    }

  // Add X2 interfaces
  mmwaveHelper->AddX2Interface (lteEnbNodes, mmWaveEnbNodes);

  // Manual attachment
    mmwaveHelper->AttachToClosestEnb (mcUeDevs, mmWaveEnbDevs, lteEnbDevs);

  // Install and start applications on UEs and remote host
    uint16_t sinkPort = 8080;

    Address sinkAddress (InetSocketAddress (ueIpIface.GetAddress (0), sinkPort));
    PacketSinkHelper packetSinkHelper ("ns3::UdpSocketFactory", InetSocketAddress (Ipv4Address::GetAny (), sinkPort));
    ApplicationContainer sinkApps = packetSinkHelper.Install (ueNodes.Get (0));
    sink = StaticCast<PacketSink> (sinkApps.Get (0));

    Ptr<Socket> ns3UdpSocket = Socket::CreateSocket (remoteHostContainer.Get (0), UdpSocketFactory::GetTypeId ());

    Ptr<MyApp> app = CreateObject<MyApp> ();
    app->Setup (ns3UdpSocket, sinkAddress, 1400, 50000000, DataRate ("1000Mb/s"));
    remoteHostContainer.Get (0)->AddApplication (app);

//    AsciiTraceHelper asciiTraceHelper;
//    Ptr<OutputStreamWrapper> stream1 = asciiTraceHelper.CreateFileStream ("mmWave-tcp-window.txt");
//    ns3TcpSocket->TraceConnectWithoutContext ("CongestionWindow", MakeBoundCallback (&CwndChange, stream1));
//
//    Ptr<OutputStreamWrapper> stream4 = asciiTraceHelper.CreateFileStream ("mmWave-tcp-rtt-newreno.txt");
//    ns3TcpSocket->TraceConnectWithoutContext ("RTT", MakeBoundCallback (&RttChange, stream4));
//
//    Ptr<OutputStreamWrapper> stream2 = asciiTraceHelper.CreateFileStream ("mmWave-tcp-data-newreno.txt");
//    sinkApps.Get (0)->TraceConnectWithoutContext ("Rx",MakeBoundCallback (&Rx, stream2));

    // ueDevs.Get (0)->TraceConnectWithoutContext ("PhyRxDrop", MakeCallback (&RxDrop));
    // Config::ConnectWithoutContext("/NodeList/*/DeviceList/*/$ns3::WifiNetDevice/Phy/PhyRxDrop ", ssMakeCallback(&PhyDrop));

    app->SetStartTime (Seconds (0.1));
    app->SetStopTime (Seconds (simTime));

    sinkApps.Start (Seconds (0.));
    sinkApps.Stop (Seconds (simTime));
    // uint16_t dlPort = 1234;
    // uint16_t ulPort = 2000;
    // ApplicationContainer clientApps;
    // ApplicationContainer serverApps;
    // bool dl = 1;
    // bool ul = 0;
    //          UdpServerHelper dlPacketSinkHelper (dlPort);
    //         dlPacketSinkHelper.SetAttribute ("PacketWindowSize", UintegerValue (256));
    //         serverApps.Add (dlPacketSinkHelper.Install (ueNodes.Get (0)));
    //         sink = StaticCast<PacketSink> (serverApps.Get (0));

    //         // Simulator::Schedule(MilliSeconds(20), &PrintLostUdpPackets, DynamicCast<UdpServer>(serverApps.Get(serverApps.GetN()-1)), lostFilename);

    //         UdpClientHelper dlClient (ueIpIface.GetAddress (0), dlPort);
    //         dlClient.SetAttribute ("Interval", TimeValue (MicroSeconds (20)));
    //         dlClient.SetAttribute ("MaxPackets", UintegerValue (0xFFFFFFFF));
    //         clientApps.Add (dlClient.Install (remoteHostContainer.Get (0)));

    // serverApps.Start (Seconds (0.));
    // clientApps.Start (Seconds (0.));
    // clientApps.Stop (Seconds (simStopTime));

    Simulator::Schedule (Seconds (0.2), &CalculateThroughput);

    Simulator::Schedule (Seconds (0.2), &ChangeSpeed, ueNodes.Get (0), Vector (ueSpeed, 0, 0)); // start UE movement after Seconds(0.5)
    // Simulator::Schedule (Seconds (1.5), &ChangeSpeed, ueNodes.Get (0), Vector (-ueSpeed, 0, 0)); // start UE movement after Seconds(0.5)
    // Simulator::Schedule (Seconds (6), &ChangeSpeed, ueNodes.Get (0), Vector (ueSpeed, 0, 0)); // start UE movement after Seconds(0.5)
    // Simulator::Schedule (Seconds (simStopTime - 1), &ChangeSpeed, ueNodes.Get (0), Vector (0, 0, 0)); // start UE movement after Seconds(0.5)
    // Simulator::Schedule (Seconds (3), &ChangeSpeed, ueNodes.Get (0), Vector (0, 0, 0));

    double numPrints = 10;
    for (int i = 0; i < numPrints; i++)
    {
        Simulator::Schedule (Seconds (i * simTime / numPrints), &PrintPosition, ueNodes.Get (0));
        // Simulator::Schedule (Seconds (i * simStopTime / numPrints), &PrintDistance, ueNodes.Get (0), allEnbNodes);
    }

     BuildingsHelper::MakeMobilityModelConsistent ();
    Config::Set ("/NodeList/*/DeviceList/*/TxQueue/MaxSize", QueueSizeValue (QueueSize ("1000000p")));
    mmwaveHelper->EnableTraces ();
    bool print = false;
    if (print)
    {
        PrintGnuplottableBuildingListToFile ("buildings.txt");
        PrintGnuplottableUeListToFile ("ues.txt");
        PrintGnuplottableEnbListToFile ("enbs.txt");
    }
    else
    {
        Simulator::Stop (Seconds (simTime));
        Simulator::Run ();
    }

    // p2ph.EnablePcapAll ("mmwave-sgi-capture");

    Simulator::Destroy ();
    myFile.close();
    return 0;
}
