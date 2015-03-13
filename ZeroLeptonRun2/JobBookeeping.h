#ifndef ZeroLeptonRun2_JobBookeeping_H_
#define ZeroLeptonRun2_JobBookeeping_H_

#include "cafe/Processor.h"

class TH1D;
class TFile;
class TObjArray;
#include <vector>
#include <map>
#include <string>

class JobBookeeping : public cafe::Processor {

public:
  JobBookeeping(const char *name);
  void begin();
  void inputFileOpened(TFile *file);
  void inputFileClosing(TFile *file);
  bool processEvent(xAOD::TEvent& event);
  void finish();

private:
  TH1D* m_counter;
  TObjArray* m_fileInfos;
  unsigned int m_eventCounter;
  std::vector<unsigned int> m_eventsPerFile;
  std::vector<std::string> m_openedFiles;
  std::vector<std::string> m_closedFiles;
  std::map<std::string, std::string> m_fileCatalog;

public:
  ClassDef(JobBookeeping,0);
};

#endif // ZeroLeptonRun2_JobBookeeping_H_
