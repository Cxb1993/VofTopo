#ifndef TIME_MEASURE_H
#define TIME_MEASURE_H

#include <string>
#include <chrono>
#include <vector>

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;
typedef std::chrono::nanoseconds Nanoseconds;
typedef std::chrono::milliseconds Milliseconds;
typedef std::chrono::seconds Seconds;

/*  #define STOPWATCH(fn,unit) {                                            \ */
/*   TimePoint start = Clock::now();                                       \ */
/*   fn;                                                                   \ */
/*   TimePoint end = Clock::now();                                         \ */
/*   unit delta = std::chrono::duration_cast<unit>(end-start);             \ */
/*   std::cout << "    " << #fn << " took " << delta.count() << " " << #unit << std::endl; \ */
/* } */

class TimeMeasure {
 public:

  int t0, t1;
  TimePoint start_all;
  TimePoint end_all;
  
  TimePoint start_seeding;
  TimePoint end_seeding;


  TimePoint start_boundary;
  TimePoint end_boundary;

  TimePoint start_components;
  TimePoint end_components;

  TimePoint start_assignment;
  TimePoint end_assignment;
  
  std::vector<TimePoint> start_oneStep;
  std::vector<TimePoint> end_oneStep;

  std::vector<TimePoint> start_interm_boundary;
  std::vector<TimePoint> end_interm_boundary;

  std::vector<TimePoint> start_advection;
  std::vector<TimePoint> end_advection;

  std::string wallClockStartTime;

  static TimeMeasure* GetInstance();

  static void StoreWallTime();

 protected:
  TimeMeasure() {}

 private:
  static TimeMeasure* instance;
};

void writeTimings(std::string filename, const TimeMeasure* timeMeasure);

#endif//TIME_MEASURE_H
