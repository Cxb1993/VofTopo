#include "time_measure.h"

#include <iostream>
#include <fstream>
#include <ctime>

TimeMeasure* TimeMeasure::instance = 0;

TimeMeasure* TimeMeasure::GetInstance() {
  if (!instance) {
    instance = new TimeMeasure();

    instance->start_oneStep.clear();
    instance->end_oneStep.clear();
    instance->start_interm_boundary.clear();
    instance->end_interm_boundary.clear();
    instance->start_advection.clear();
    instance->end_advection.clear();
  }
  return instance;
}

void TimeMeasure::StoreWallTime() {
  time_t rawtime;
  struct tm * timeinfo;
  char buffer [80];
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,80,"%F_%H:%M:%S",timeinfo);
  
  instance->wallClockStartTime = buffer;
}

void writeTimings(std::string filename, const TimeMeasure* tm)
{
  std::ofstream file(filename, std::ofstream::out);

  file << "timesteps," << tm->t0 << "," << tm->t1 << std::endl;

  Milliseconds dur_all_ms =
    std::chrono::duration_cast<Milliseconds>(tm->end_all-tm->start_all);
  file << "all," << dur_all_ms.count() << std::endl;

  Milliseconds dur_boundary_ms =
    std::chrono::duration_cast<Milliseconds>(tm->end_boundary-tm->start_boundary);
  file << "boundary," << dur_boundary_ms.count() << std::endl;

  Milliseconds dur_components_ms =
    std::chrono::duration_cast<Milliseconds>(tm->end_components-tm->start_components);
  file << "components," << dur_components_ms.count() << std::endl;

  Milliseconds dur_assignment_ms =
    std::chrono::duration_cast<Milliseconds>(tm->end_assignment-tm->start_assignment);
  file << "assignment," << dur_assignment_ms.count() << std::endl;

  if (tm->start_oneStep.size() == tm->end_oneStep.size()) {
    for (int i = 0; i < tm->start_oneStep.size(); ++i) {
      Milliseconds dur_step_ms =
	std::chrono::duration_cast<Milliseconds>(tm->end_oneStep[i]-tm->start_oneStep[i]);
      file << "oneStep," << i << "," << dur_step_ms.count() << std::endl;
    }
  }
  if (tm->start_advection.size() == tm->end_advection.size()) {
    for (int i = 0; i < tm->start_advection.size(); ++i) {
      Milliseconds dur_step_ms =
	std::chrono::duration_cast<Milliseconds>(tm->end_advection[i]-tm->start_advection[i]);
      file << "advection," << i << "," << dur_step_ms.count() << std::endl;
    }
  }
  if (tm->start_interm_boundary.size() == tm->end_interm_boundary.size()) {
    for (int i = 0; i < tm->start_interm_boundary.size(); ++i) {
      Milliseconds dur_step_ms =
	std::chrono::duration_cast<Milliseconds>(tm->end_interm_boundary[i]-tm->start_interm_boundary[i]);
      file << "interm_boundary," << i << "," << dur_step_ms.count() << std::endl;
    }
  }

  file.close();
}
