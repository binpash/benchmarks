FROM hs-riker

RUN apt-get update && apt-get install -y git clang cmake

WORKDIR /

RUN git clone https://github.com/llvm/llvm-project.git
WORKDIR /llvm-project
RUN git checkout d28af7c654d8db0b68c175db5ce212d74fb5e9bc
RUN rm -rf llvm-build && mkdir llvm-build

WORKDIR /llvm-project/llvm-build
RUN cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_ENABLE_BINDINGS=0 ../llvm

COPY llvm_build.sh /llvm-project/llvm-build/llvm_build.sh
COPY llvm_build.sh /llvm-project/llvm-build/Rikerfile
RUN chmod +x /llvm-project/llvm-build/Rikerfile
RUN chmod +x /llvm-project/llvm-build/llvm_build.sh