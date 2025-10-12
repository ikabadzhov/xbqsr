cmake_minimum_required(VERSION 3.22)
SET(CMAKE_C_COMPILER icx)
SET(CMAKE_CXX_COMPILER icpx)

project(bqsr)

SET(CMAKE_CXX_STANDARD 20)


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -g -O2")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsycl -fsycl-targets=nvptx64-nvidia-cuda,spir64_x86_64 -Xsycl-target-backend=nvptx64-nvidia-cuda --cuda-gpu-arch=sm_75 -Wall")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsycl -fsycl-targets=nvptx64-nvidia-cuda,spir64_x86_64 -Xsycl-target-backend=nvptx64-nvidia-cuda --cuda-gpu-arch=sm_89 -Wall")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lOpenCL -lsycl")


find_package(ZLIB REQUIRED)
find_library(LIBDEFLATE_LIBRARY deflate REQUIRED)
find_package(Boost 1.63.0 COMPONENTS program_options REQUIRED)
find_package(TBB)
include_directories( ${Boost_INCLUDE_DIR})

add_executable(bqsr main.cpp
		sam_data_structures/sam_fields/Sequence.cpp
		input_output/utils/StringScanner.cpp
		sam_data_structures/SamHeader.cpp
		input_output/alignment_io/AlignmentIO.cpp
        Utils.cpp
		sam_data_structures/sam_fields/Cigar.cpp file_format/Global.cpp filters/Filter.cpp
		filters/FilterAlignment.cpp Utils.cpp
		bqsr/BQSR.cpp bqsr/BQSR.h intervals/Interval.cpp intervals/Interval.h
        filters/UtilsFilters.cpp filters/UtilsFilters.h
		file_format/MappedFasta.cpp file_format/MappedFasta.h
		bqsr/TablesBQSR.h bqsr/TablesBQSR.cpp
		sam_data_structures/SamAlignments.cpp sam_data_structures/SamAlignments.h
		sam_data_structures/sam_fields/OptionalFields.cpp sam_data_structures/sam_fields/OptionalFields.h
		sam_data_structures/sam_fields/TmpFields.cpp sam_data_structures/sam_fields/TmpFields.h
		data_types/CompressedStringType.cpp data_types/CompressedStringCollection.h
		data_types/StringCollection.cpp data_types/StringCollection.h
		data_types/VectorCollection.cpp data_types/VectorCollection.h
		data_types/TagsCollection.h data_types/BitVector.cpp data_types/BitVector.h
		data_types/TagsCollection.cpp input_output/sam_io/SamWriter.cpp
		data_types/LIBIDStringMapper.cpp
		data_types/RGStringMapper.cpp
		input_output/sam_io/SamWriter.h input_output/sam_io/SamReader.cpp
		input_output/sam_io/SamReader.h input_output/bam_io/BamReader.cpp
		input_output/bam_io/BamReader.h input_output/bam_io/BamWriter.cpp
		input_output/bam_io/BamWriter.h
		bam_data_structures/BamHeader.cpp
		bam_data_structures/BamHeader.h
		input_output/alignment_io/AlignmentReaderFactory.cpp
		input_output/alignment_io/AlignmentReaderFactory.h
        input_output/streams/CompressedInputStream.cpp
        input_output/streams/CompressedInputStream.h
		input_output/streams/PlainInputStream.cpp
		input_output/streams/PlainInputStream.h
        input_output/streams/SamStream.cpp
        input_output/streams/SamStream.h
		input_output/streams/InputStream.h
        input_output/streams/BufferReader.cpp
        input_output/streams/BufferReader.h
		input_output/streams/BufferWriter.cpp
		input_output/streams/BufferWriter.h
		input_output/utils/ByteScanner.cpp
		data_types/LinkedHashMap.h
		sam_data_structures/sam_header_fields/ReadGroupField.cpp
		sam_data_structures/sam_header_fields/ReadGroupField.h
		sam_data_structures/sam_header_fields/SequenceField.h
		sam_data_structures/sam_header_fields/ProgramField.h
		input_output/parser/ParsedHeaderLine.cpp
		input_output/parser/ParsedHeaderLine.h
		input_output/parser/SamHeaderParser.cpp
		input_output/parser/SamHeaderParser.h
		sam_data_structures/sam_header_fields/SamHeaderRecord.h
		sam_data_structures/sam_header_fields/FileMetadata.h
		sam_data_structures/sam_header_fields/SequenceField.cpp
		sam_data_structures/sam_header_fields/ProgramField.cpp
		data_types/QnameStringMapper.cpp
		input_output/parser/SamAlignmentParser.cpp
		input_output/parser/SamAlignmentParser.h
		input_output/formatter/SamHeaderFormatter.cpp
		input_output/formatter/SamHeaderFormatter.h
		input_output/formatter/SamAlignmentFormatter.cpp
		input_output/formatter/SamAlignmentFormatter.h
		sam_data_structures/sam_header_fields/SamHeaderRecord.cpp
		input_output/parser/BamAlignmentParser.cpp
		input_output/parser/BamAlignmentParser.h
		input_output/formatter/BamAlignmentFormatter.cpp
		input_output/formatter/BamAlignmentFormatter.h
		input_output/formatter/BamHeaderFormatter.cpp
		input_output/formatter/BamHeaderFormatter.h
		input_output/alignment_io/AlignmentWriterFactory.cpp
		input_output/alignment_io/AlignmentWriterFactory.h
		input_output/streams/OutputStream.h
		input_output/streams/PlainOutputStream.cpp
		input_output/streams/PlainOutputStream.h
		input_output/streams/CompressedOutputStream.cpp
		input_output/streams/CompressedOutputStream.h
		input_output/utils/ByteFormatter.cpp
		input_output/utils/ByteFormatter.h
		data_types/SparseVector.cpp
		data_types/SparseVector.h
		sam_data_structures/sam_fields/Tag.h
		sam_data_structures/sam_fields/Tag.cpp
		input_output/utils/Constant.h
		timer/Timer.cpp
		sycl_parallel_library/device_utils/types.cpp
		sycl_parallel_library/device_utils/meta.h
		sycl_parallel_library/device_utils/operator.h
		sycl_parallel_library/device_utils/types.h
		sycl_parallel_library/merge/merge.cpp
		sycl_parallel_library/merge/merge.h
		sycl_parallel_library/merge/external_merge.h
		sycl_parallel_library/merge/search.h
		sycl_parallel_library/sorting/mergesort.h
		sycl_parallel_library/sorting/cta_merge.h
		sycl_parallel_library/sorting/cta_mergesort.cpp
		sycl_parallel_library/sorting/cta_mergesort.h
		sycl_parallel_library/sorting/sort_network.h
		sycl_parallel_library/memory/context.h
		sycl_parallel_library/memory/loadstore.h
		sycl_parallel_library/memory/memory.h
		sycl_parallel_library/load_balancing/load_balance.h
		sycl_parallel_library/intervalsmove/intervalmove.h
		bqsr/compute_bqsr_tables/ComputeBQSRTables.cpp
		bqsr/compute_bqsr_tables/ComputeBQSRTables.h
		bqsr/compute_bqsr_tables/ComputeBQSRTablesKernelData.h
		bqsr/compute_vectors_for_covariates/ComputeCovariateVectorsFunctions.cpp
		bqsr/compute_vectors_for_covariates/ComputeCovariateVectorsFunctions.h
		bqsr/compute_vectors_for_covariates/ComputeCovariateVectorsData.h
		bqsr/compute_vectors_for_covariates/ComputeCovariateVectorsData.cpp
		bqsr/UtilsBQSR.cpp
		input_output/parser/IOBatch.cpp
		input_output/parser/IOBatch.h
		Config.cpp
		Config.h
		sycl_parallel_library/segmented_scan/segscan.h
		bqsr/apply_bqsr/ApplyBqsrKernelData.cpp
		bqsr/apply_bqsr/ApplyBqsrKernelData.h
		bqsr/apply_bqsr/ApplyBqsrFunctions.cpp
		bqsr/apply_bqsr/ApplyBqsrFunctions.h
		tileinfo/tileInfo.cpp
		tileinfo/tileInfo.h
		input_output/bam_io/BamUtilsIO.h
		timer/Profiler.cpp
		timer/Profiler.h
		timer/ProfilerIO.h
)

target_link_libraries(bqsr TBB::tbb /opt/intel/oneapi/advisor/latest/lib64/libittnotify.a pthread dl atomic ${Boost_LIBRARIES} ZLIB::ZLIB ${LIBDEFLATE_LIBRARY})