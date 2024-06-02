/* ==========================================================================
 *  Copyright (C) KU Leuven, 2019
 * 
 *  This file is licensed under the terms of the 3-Clause BSD License. You may
 *  not use this file except in compliance with the 3-Clause BSD License. You
 *  should have received a copy of the 3-clause BSD license along with this
 *  file. If not, see <https://opensource.org/licenses/BSD-3-Clause>.
 * ==========================================================================*/

#ifndef TICTOC_H
#define TICTOC_H

#include <iostream>
#include <chrono>

/* TICTOC is a simple timer based on std::chrono::steady_clock
 * 
 * Author(s): F.P.B. (KU Leuven)
 * 
 * Mode 1:
 *   TicToc::tic();
 *   ...
 *   TicToc::toc();
 * 
 * Mode 2:
 *   t_point t0 = TicToc::tic(0);
 *   ...
 *   TicToc::toc(t0, "optional_output_text");
 */

typedef std::chrono::steady_clock::time_point t_point;
static t_point t1 = std::chrono::steady_clock::now();

namespace TicToc {
  // Declarations
  inline void tic ();
  inline t_point tic (bool dummy);

  inline void toc ();
  inline void toc (std::string outputText);
  inline void toc (t_point t0);
  inline void toc (t_point t0, std::string outputText);

  inline void elapsedTime (std::chrono::duration<double> t, std::string outputText);

  // Definitions
  void tic () { t1 = std::chrono::steady_clock::now(); }
  t_point tic (bool dummy) { return std::chrono::steady_clock::now(); }

  void toc () {
    t_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_span = 
      std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    elapsedTime(time_span, "");
  }

  void toc (std::string outputText) {
    t_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_span = 
      std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    elapsedTime(time_span, outputText);
  }

  void toc (t_point t0) {
    t_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_span = 
      std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t0);
    elapsedTime(time_span, "");
  }
  void toc (t_point t0, std::string outputText) {
    t_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_span = 
      std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t0);
    elapsedTime(time_span, outputText);
  }
  
  void elapsedTime (std::chrono::duration<double> t, std::string outputText) {
    if (outputText.empty())
      outputText = "  Elapsed time: ";
    
    if ( t.count() <= 60 )
      std::cout << outputText << t.count() << " seconds" << std::endl;
    else if ( t.count() > 60 && t.count() <= 3600 )
      std::cout << outputText << t.count()/60 << " minutes" << std::endl;
    else
      std::cout << outputText << t.count()/3600 << " hours" << std::endl;
  }

}
#endif