class Spdlog < Formula
  desc "Super fast C++ logging library"
  homepage "https://github.com/gabime/spdlog"
  url "https://github.com/gabime/spdlog/archive/v1.8.5.tar.gz"
  sha256 "944d0bd7c763ac721398dca2bb0f3b5ed16f67cef36810ede5061f35a543b4b8"
  license "MIT"
  head "https://github.com/gabime/spdlog.git", branch: "v1.x"

  depends_on "cmake" => :build
  depends_on "fmt"

  def install
    ENV.cxx11

    inreplace "include/spdlog/tweakme.h", "// #define SPDLOG_FMT_EXTERNAL", "#define SPDLOG_FMT_EXTERNAL"

    mkdir "spdlog-build" do
      args = std_cmake_args + %W[
        -Dpkg_config_libdir=#{lib}
        -DSPDLOG_BUILD_BENCH=OFF
        -DSPDLOG_BUILD_TESTS=OFF
        -DSPDLOG_FMT_EXTERNAL=ON
      ]
      system "cmake", "..", "-DSPDLOG_BUILD_SHARED=ON", *args
      system "make", "install"
      system "make", "clean"
      system "cmake", "..", "-DSPDLOG_BUILD_SHARED=OFF", *args
      system "make"
      lib.install "libspdlog.a"
    end
  end

  test do
    (testpath/"test.cpp").write <<~EOS
      #include "spdlog/sinks/basic_file_sink.h"
      #include <iostream>
      #include <memory>
      int main()
      {
        try {
          auto console = spdlog::basic_logger_mt("basic_logger", "#{testpath}/basic-log.txt");
          console->info("Test");
        }
        catch (const spdlog::spdlog_ex &ex)
        {
          std::cout << "Log init failed: " << ex.what() << std::endl;
          return 1;
        }
      }
    EOS

    system ENV.cxx, "-std=c++11", "test.cpp", "-I#{include}", "-L#{Formula["fmt"].opt_lib}", "-lfmt", "-o", "test"
    system "./test"
    assert_predicate testpath/"basic-log.txt", :exist?
    assert_match "Test", (testpath/"basic-log.txt").read
  end
end