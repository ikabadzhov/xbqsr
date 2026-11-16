# X-BQSR: Holistic Acceleration of Base Quality Recalibration for Scalable Genomic Analysis

X-BQSR is a high-throughput, cross-vendor GPU-accelerated implementation of the GATK Base Quality Score Recalibration (BQSR) workflow using SYCL and [AdaptiveCpp](https://github.com/AdaptiveCpp/AdaptiveCpp). It corrects systematic sequencing errors by analyzing quality score patterns across read groups, cycle positions, and dinucleotide contexts, generating recalibrated BAM files with bitwise and statistical parity against GATK.

---

## Key Features

- **Cross-Vendor GPU Acceleration**: Native high-performance execution on both **NVIDIA** (CUDA / NVPTX) and **AMD** (ROCm / HIP / AMDGPU) GPUs.
- **AdaptiveCpp Multi-Backend Compilation**: Single-Pass Compiler (SSCP) `generic` target enables a single binary to JIT-compile and execute across heterogeneous GPU architectures without re-compilation.
- **Two-Stage Workflow**:
  1. `recalibrate` – Ingests aligned BAM, reference FASTA, and known sites VCFs to compute empirical quality score tables.
  2. `apply` – Applies the recalibration tables to adjust base quality scores in the BAM stream.
- **nf-core / Sarek Ready**: Seamlessly integrated into nf-core/sarek as an accelerated BQSR engine (`--bqsr_tool xbqsr`).
- **Hardware Requirement**: X-BQSR is designed strictly for GPU accelerators. CPU execution is not supported.

---

## Hardware Requirements

- **GPU Accelerators**:
  - **NVIDIA**: Compute Capability &ge; `sm_70` (Volta, Turing, Ampere, Ada Lovelace, Hopper). Requires NVIDIA driver and CUDA Toolkit.
  - **AMD**: CDNA / RDNA architectures (e.g. Instinct MI100/MI200/MI300, Radeon Pro). Requires ROCm &ge; 5.7.
- **CPU Fallback**: X-BQSR kernels rely on GPU Cooperative Thread Array (CTA) work-group primitives. CPU execution is **not supported**; the binary explicitly requires a GPU accelerator and will exit with code 1 if no compatible GPU is detected. For CPU-only clusters, use GATK (`--bqsr_tool gatk`).

---

## Building from Source

### Prerequisites

- **C++ Compiler**: Clang &ge; 15 or GCC &ge; 11
- **AdaptiveCpp**: &ge; 24.02.0 installed with CUDA and/or ROCm backends enabled
- **Libraries**:
  - Boost (`libboost-program-options-dev`, `libboost-fiber-dev`, `libboost-context-dev`, `libboost-test-dev`)
  - Intel TBB (`libtbb-dev`)
  - libdeflate (`libdeflate-dev`)
  - zlib (`zlib1g-dev`)
- **Build Tools**: CMake &ge; 3.20, Ninja

### Native Build with AdaptiveCpp (Default)

```bash
# Clone the repository
git clone https://gitlab.eurecom.fr/dislab/bqsr.git
cd bqsr

# Configure using AdaptiveCpp (generic SSCP target by default)
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=acpp \
  -DBQSR_USE_ADAPTIVECPP=ON \
  -DACPP_TARGETS=generic

# Build
cmake --build build -j$(nproc)
```

The resulting executable will be generated at `build/bqsr`.

### Explicit Ahead-Of-Time (AOT) Targets (Optional)

To compile for specific GPU architectures instead of using the generic JIT:

```bash
# For NVIDIA L40S / Ada Lovelace:
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=acpp -DACPP_TARGETS="cuda:sm_89"

# For AMD Instinct MI210:
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=acpp -DACPP_TARGETS="hip:gfx90a"

# Combined NVIDIA + AMD binary:
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=acpp -DACPP_TARGETS="cuda:sm_80,sm_89;hip:gfx90a"
```

---

## Usage

### Check Version & Help
```bash
./build/bqsr -V
./build/bqsr --help
./build/bqsr recalibrate --help
./build/bqsr apply --help
```

### 1. Recalibrate Stage
Generate recalibration tables from aligned reads:
```bash
./build/bqsr recalibrate \
  --bam-in sample.bam \
  --recal-out sample.recal.table \
  --ref reference.fasta \
  --known-sites dbsnp.vcf.gz,mills_1000G.vcf.gz
```

*Note: The input BAM file must have an accompanying index (`.bam.bai` or `.bai`) in the same directory.*

### 2. Apply Stage
Apply the recalibration table to adjust base qualities:
```bash
./build/bqsr apply \
  --bam-in sample.bam \
  --bam-out sample.recal.bam \
  --recal sample.recal.table
```

### Command-Line Arguments

| Option | Command | Description |
|---|---|---|
| `--bam-in <file>` | Both | Input BAM file with aligned reads (requires `.bai` index) |
| `--recal-out <file>` | `recalibrate` | Output recalibration table file |
| `--ref <file>` | `recalibrate` | Reference genome FASTA file |
| `--known-sites <files>` | `recalibrate` | Comma-delimited list of known variant VCF/VCF.GZ files |
| `--bam-out <file>` | `apply` | Output recalibrated BAM file |
| `--recal <file>` | `apply` | Recalibration table file from the `recalibrate` stage |
| `--gpus <num\|all>` | Both | GPU selection: specific index (e.g. `0`), list (`0,1`), or `all` |
| `--verbose, -v` | Both | Enable verbose logging |
| `--version, -V` | Global | Print version information |
| `--help, -h` | Global | Display help message |

---

## Container Usage

### Docker (NVIDIA)
```bash
docker run --rm --gpus all \
  -v /path/to/data:/data \
  gitlab.eurecom.fr:5050/dislab/bqsr:feature-sarek-integration \
  recalibrate --bam-in /data/input.bam --recal-out /data/recal.table \
              --ref /data/ref.fasta --known-sites /data/known.vcf.gz
```

### Apptainer / Singularity
```bash
# NVIDIA:
apptainer exec --nv docker://gitlab.eurecom.fr:5050/dislab/bqsr:feature-sarek-integration \
  bqsr recalibrate --bam-in input.bam --recal-out recal.table --ref ref.fasta --known-sites sites.vcf.gz

# AMD ROCm:
apptainer exec --rocm docker://gitlab.eurecom.fr:5050/dislab/bqsr:feature-sarek-integration \
  bqsr recalibrate --bam-in input.bam --recal-out recal.table --ref ref.fasta --known-sites sites.vcf.gz
```

---

## nf-core / Sarek Integration

X-BQSR is integrated into nf-core/sarek as an accelerated base quality recalibration alternative:

```bash
nextflow run nf-core/sarek \
  -profile docker,gpu \
  --input samplesheet.csv \
  --outdir results/ \
  --bqsr_tool xbqsr
```

---

## Citation

If you use X-BQSR in your research, please cite:

```bibtex
@InProceedings{Eugenio2027xbqsr,
  author    = "Marinelli, Eugenio and Kabadzhov, Ivan Donchev and Appuswamy, Raja",
  editor    = "Torquati, Massimo and Mencagli, Gabriele and Cardellini, Valeria and Antelmi, Alessia and Medi{\'{c}}, Doriana",
  title     = "X-BQSR: Holistic Acceleration of Base Quality Recalibration for Scalable Genomic Analysis",
  booktitle = "Euro-Par 2026: Parallel Processing",
  year      = "2027",
  publisher = "Springer Nature Switzerland",
  address   = "Cham",
  pages     = "523--538",
  abstract  = "With the continued decline in sequencing costs, genomics-based precision medicine is gaining increasing importance. The Genome Analysis Toolkit (GATK) is the standard framework for converting raw genomic data into actionable variants, but one of its most expensive steps is base quality score recalibration (BQSR), which corrects systematic sequencing errors. Despite significant work on accelerating other parts of the variant calling pipeline, BQSR has not been systematically explored on GPUs.",
  isbn      = "978-3-032-35248-4"
}
```

---

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details.
