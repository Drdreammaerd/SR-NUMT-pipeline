#!/bin/bash
# Summarize NUMT pipeline benchmark results from LSF job logs
# Usage: bash summarize_benchmarks.sh /path/to/logs_dir

LOG_DIR="${1:-/Volumes/jin810/Active/testing/yung-chun/AI-develop/SMAHT_DONOR_NUMT/logs}"

echo "============================================================"
echo " NUMT Pipeline — Benchmark Summary"
echo "============================================================"

# --- Orchestrator ---
echo ""
echo ">>> ORCHESTRATOR (Total Wall-Clock)"
grep "Started at\|Terminated at" "${LOG_DIR}/orchestrator.log" | tail -2
ORCH_RT=$(grep "Run time" "${LOG_DIR}/orchestrator.log" | tail -1 | sed 's/[^0-9.]//g')
echo "    Total Wall-Clock: ${ORCH_RT} sec ($(echo "scale=1; ${ORCH_RT}/3600" | bc) hours)"

# --- Per-rule summary ---
echo ""
echo ">>> TOP 15 SLOWEST JOBS (by Run time)"
echo "Run_Time(s)  CPU_Time(s)  Max_Mem(MB)  Job_Name"
echo "----------  -----------  -----------  --------"

for f in "${LOG_DIR}"/run_dinumt_split_*.out "${LOG_DIR}"/split_bam_*.out "${LOG_DIR}"/extract_insert_size_*.out "${LOG_DIR}"/merge_vcfs_*.out "${LOG_DIR}"/generate_manifest_*.out "${LOG_DIR}"/build_final_report_*.out; do
    [ -f "$f" ] || continue
    rt=$(grep "Run time" "$f" 2>/dev/null | sed 's/[^0-9.]//g')
    ct=$(grep "CPU time" "$f" 2>/dev/null | sed 's/[^0-9.]//g')
    mm=$(grep "Max Memory" "$f" 2>/dev/null | sed 's/[^0-9.]//g' | head -1)
    name=$(basename "$f" .out)
    [ -n "$rt" ] && echo "$rt  $ct  $mm  $name"
done | sort -rn | head -15

# --- Aggregate by rule ---
echo ""
echo ">>> AGGREGATE BY RULE (count, avg wall-time, max wall-time)"

for rule in run_dinumt_split split_bam extract_insert_size merge_vcfs generate_manifest build_final_report; do
    count=0; total=0; maxrt=0
    for f in "${LOG_DIR}"/${rule}_*.out; do
        [ -f "$f" ] || continue
        rt=$(grep "Run time" "$f" 2>/dev/null | sed 's/[^0-9.]//g')
        [ -z "$rt" ] && continue
        count=$((count + 1))
        total=$(echo "$total + $rt" | bc)
        if [ "$(echo "$rt > $maxrt" | bc)" -eq 1 ]; then maxrt=$rt; fi
    done
    [ "$count" -gt 0 ] && avg=$(echo "scale=1; $total / $count" | bc) || avg=0
    echo "  ${rule}: ${count} jobs, avg=${avg}s, max=${maxrt}s, total_compute=${total}s ($(echo "scale=1; $total/3600" | bc)h)"
done

echo ""
echo "============================================================"
