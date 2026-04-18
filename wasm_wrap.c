#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <stdarg.h>
#include <limits.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include <emscripten/emscripten.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ketopt.h"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
    struct rlimit r;
    getrlimit(RLIMIT_AS, &r);
    r.rlim_cur = r.rlim_max;
    setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

static ko_longopt_t long_options[] = {
    { "bucket-bits",    ko_required_argument, 300 },
    { "mb-size",        ko_required_argument, 'K' },
    { "seed",           ko_required_argument, 302 },
    { "no-kalloc",      ko_no_argument,       303 },
    { "print-qname",    ko_no_argument,       304 },
    { "no-self",        ko_no_argument,       'D' },
    { "print-seeds",    ko_no_argument,       306 },
    { "max-chain-skip", ko_required_argument, 307 },
    { "min-dp-len",     ko_required_argument, 308 },
    { "print-aln-seq",  ko_no_argument,       309 },
    { "splice",         ko_no_argument,       310 },
    { "cost-non-gt-ag", ko_required_argument, 'C' },
    { "no-long-join",   ko_no_argument,       312 },
    { "sr",             ko_optional_argument, 313 },
    { "frag",           ko_required_argument, 314 },
    { "secondary",      ko_required_argument, 315 },
    { "cs",             ko_optional_argument, 316 },
    { "end-bonus",      ko_required_argument, 317 },
    { "no-pairing",     ko_no_argument,       318 },
    { "splice-flank",   ko_required_argument, 319 },
    { "idx-no-seq",     ko_no_argument,       320 },
    { "end-seed-pen",   ko_required_argument, 321 },
    { "for-only",       ko_no_argument,       322 },
    { "rev-only",       ko_no_argument,       323 },
    { "heap-sort",      ko_required_argument, 324 },
    { "all-chain",      ko_no_argument,       'P' },
    { "dual",           ko_required_argument, 326 },
    { "max-clip-ratio", ko_required_argument, 327 },
    { "min-occ-floor",  ko_required_argument, 328 },
    { "MD",             ko_no_argument,       329 },
    { "lj-min-ratio",   ko_required_argument, 330 },
    { "score-N",        ko_required_argument, 331 },
    { "eqx",            ko_no_argument,       332 },
    { "paf-no-hit",     ko_no_argument,       333 },
    { "split-prefix",   ko_required_argument, 334 },
    { "no-end-flt",     ko_no_argument,       335 },
    { "hard-mask-level",ko_no_argument,       336 },
    { "cap-sw-mem",     ko_required_argument, 337 },
    { "max-qlen",       ko_required_argument, 338 },
    { "max-chain-iter", ko_required_argument, 339 },
    { "junc-bed",       ko_required_argument, 340 },
    { "junc-bonus",     ko_required_argument, 341 },
    { "sam-hit-only",   ko_no_argument,       342 },
    { "chain-gap-scale",ko_required_argument, 343 },
    { "alt",            ko_required_argument, 344 },
    { "alt-drop",       ko_required_argument, 345 },
    { "mask-len",       ko_required_argument, 346 },
    { "rmq",            ko_optional_argument, 347 },
    { "qstrand",        ko_no_argument,       348 },
    { "cap-kalloc",     ko_required_argument, 349 },
    { "q-occ-frac",     ko_required_argument, 350 },
    { "chain-skip-scale",ko_required_argument,351 },
    { "print-chains",   ko_no_argument,       352 },
    { "no-hash-name",   ko_no_argument,       353 },
    { "secondary-seq",  ko_no_argument,       354 },
    { "ds",             ko_no_argument,       355 },
    { "rmq-inner",      ko_required_argument, 356 },
    { "spsc",           ko_required_argument, 357 },
    { "junc-pen",       ko_required_argument, 358 },
    { "pairing",        ko_required_argument, 359 },
    { "jump-min-match", ko_required_argument, 360 },
    { "write-junc",     ko_no_argument,       361 },
    { "pass1",          ko_required_argument, 362 },
    { "spsc-scale",     ko_required_argument, 363 },
    { "spsc0",          ko_required_argument, 364 },
    { "dbg-seed-occ",   ko_no_argument,       501 },
    { "help",           ko_no_argument,       'h' },
    { "max-intron-len", ko_required_argument, 'G' },
    { "version",        ko_no_argument,       'V' },
    { "min-count",      ko_required_argument, 'n' },
    { "min-chain-score",ko_required_argument, 'm' },
    { "mask-level",     ko_required_argument, 'M' },
    { "min-dp-score",   ko_required_argument, 's' },
    { "sam",            ko_no_argument,       'a' },
    { 0, 0, 0 }
};

static inline int64_t mm_parse_num2(const char *str, char **q)
{
    double x;
    char *p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
    else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
    else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
    if (q) *q = p;
    return (int64_t)(x + .499);
}

static inline int64_t mm_parse_num(const char *str)
{
    return mm_parse_num2(str, 0);
}

static inline void yes_or_no(mm_mapopt_t *opt, int64_t flag, int long_idx, const char *arg, int yes_to_set)
{
    if (yes_to_set) {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
        else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
        else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
    } else {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
        else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
        else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
    }
}

int minimap2_main(int argc, char *argv[]);

void print_minimap2_options(const mm_idxopt_t *ipt, const mm_mapopt_t *opt) {
    // fprintf(stdout, "\n--- Final Minimap2 Parameters ---\n");
    // fprintf(stdout, "[Indexing Options]\n");
    // fprintf(stdout, "  k-mer size (k): %d\n", ipt->k);
    // fprintf(stdout, "  Window size (w): %d\n", ipt->w);
    // fprintf(stdout, "  Homopolymer-compressed (HPC): %s\n", (ipt->flag & MM_I_HPC) ? "yes" : "no");
    
    // fprintf(stdout, "[Mapping Options]\n");
    // fprintf(stdout, "  Max gap (g): %d\n", opt->max_gap);
    // fprintf(stdout, "  Min chain count (n): %d\n", opt->min_cnt);
    // fprintf(stdout, "  Min chain score (m): %d\n", opt->min_chain_score);
    // fprintf(stdout, "  Preset is Short Read (SR): %s\n", (opt->flag & MM_F_SR) ? "yes" : "no");

    // fprintf(stdout, "[Alignment Options]\n");
    // fprintf(stdout, "  Match score (A): %d\n", opt->a);
    // fprintf(stdout, "  Mismatch penalty (B): %d\n", opt->b);
    // fprintf(stdout, "  Gap open penalty (O1,O2): %d, %d\n", opt->q, opt->q2);
    // fprintf(stdout, "  Gap extension penalty (E1,E2): %d, %d\n", opt->e, opt->e2);

    // fprintf(stdout, "[Output Options]\n");
    // fprintf(stdout, "  Output SAM format: %s\n", (opt->flag & MM_F_OUT_SAM) ? "yes" : "no");
    // fprintf(stdout, "  Output CG tag: %s\n", (opt->flag & MM_F_OUT_CG) ? "yes" : "no");
    // fprintf(stdout, "  Output CS tag: %s\n", (opt->flag & MM_F_OUT_CS) ? "yes" : "no");
    // fprintf(stdout, "-----------------------------------\n\n");
}


// 可调参数
#define MAX_SEGMENTS 65536
#define SAM_LINE_BUF 65536
#define MAX_FIELDS 64

typedef struct {
    uint32_t ref_start;
    uint32_t ref_end;
    uint32_t query_start;
    uint32_t query_end;
    int8_t strand;   // +1 or -1
    uint8_t score;   // mapq or similar
} segment_t;

typedef struct {
    uint32_t ref_length;    // 用于 x 轴归一化 (取解析到的最大 ref_end)
    uint32_t query_length;  // 用于 y 轴归一化 (取解析到的最大 query_end)
    segment_t segments[MAX_SEGMENTS];
    size_t n_segments;
} dotplot_t;

// append helper
static bool append_segment(dotplot_t *plot,
                           uint32_t ref_start, uint32_t ref_end,
                           uint32_t query_start, uint32_t query_end,
                           int8_t strand, uint8_t score) {
    if (!plot) return false;
    if (plot->n_segments >= MAX_SEGMENTS) return false;
    plot->segments[plot->n_segments++] = (segment_t){
        .ref_start = ref_start,
        .ref_end = ref_end,
        .query_start = query_start,
        .query_end = query_end,
        .strand = strand,
        .score = score
    };
    return true;
}

// 解析 CIGAR，把每个连续 M/=/X 区块作为一个 segment
// ref_pos and query_pos 是 IN/OUT，表示当前 reference/query 0-based 位置
static void advance_cigar_blocks(const char *cigar,
                                 uint32_t *ref_pos,
                                 uint32_t *query_pos,
                                 dotplot_t *plot,
                                 int8_t strand,
                                 uint8_t score) {
    if (!cigar || !ref_pos || !query_pos || !plot) return;
    uint32_t num = 0;
    for (const char *p = cigar; *p; ++p) {
        if (isdigit((unsigned char)*p)) {
            num = num * 10 + (uint32_t)(*p - '0');
            continue;
        }
        uint32_t len = (num == 0) ? 1 : num;
        switch (*p) {
            case 'M': case '=': case 'X': {
                uint32_t ref_start = *ref_pos;
                uint32_t query_start = *query_pos;
                *ref_pos += len;
                *query_pos += len;
                // append segment for this matching block
                append_segment(plot, ref_start, *ref_pos, query_start, *query_pos, strand, score);
                break;
            }
            case 'I': case 'S':
                *query_pos += len;
                break;
            case 'D': case 'N':
                *ref_pos += len;
                break;
            case 'H': case 'P':
                // H and P do not consume sequence (H does not consume query), skip
                break;
            default:
                // unknown op: ignore
                break;
        }
        num = 0;
    }
}

// 把 dotplot 写成 JSON 文件 path（覆盖写入）
static bool dotplot_to_json_file(const dotplot_t *plot, const char *out_path) {
    if (!plot || !out_path) return false;
    FILE *f = fopen(out_path, "wb");
    if (!f) return false;
    fprintf(f, "{\"ref_length\":%u,\"query_length\":%u,\"segments\":[",
            plot->ref_length, plot->query_length);
    for (size_t i = 0; i < plot->n_segments; ++i) {
        const segment_t *s = &plot->segments[i];
        fprintf(f,
            "{\"ref_start\":%u,\"ref_end\":%u,"
            "\"query_start\":%u,\"query_end\":%u,"
            "\"strand\":%d,\"score\":%u}%s",
            s->ref_start, s->ref_end,
            s->query_start, s->query_end,
            (int)s->strand, (unsigned)s->score,
            (i + 1 == plot->n_segments) ? "" : ",");
    }
    fprintf(f, "]}");
    fclose(f);
    return true;
}

// 解析一行 SAM alignment，返回 true 表示成功解析并加入 segments（否则 false）
// line: 包含末尾 '\n' 的缓冲区（不可修改之外的副作用）
// plot: dotplot 输出容器（会累积 segments）
// track_max: 输出时更新最大 ref/query 端点（用于计算长度）
static bool parse_sam_alignment_line(char *line, dotplot_t *plot,
                                     uint32_t *out_max_ref, uint32_t *out_max_query) {
    if (!line || !plot) return false;
    // Split by tabs to get fields. We only need fields 2(FLAG),3(RNAME),4(POS),5(MAPQ),6(CIGAR),10(SEQ)
    // Use in-place tokenization (will modify line)
    char *fields[12];
    int fcount = 0;
    char *p = line;
    char *start = p;
    while (*p && fcount < 11) {
        if (*p == '\t' || *p == '\n' || *p == '\r') {
            *p = '\0';
            fields[fcount++] = start;
            p++;
            start = p;
        } else {
            p++;
        }
    }
    // capture last field if not counted
    if (fcount < 11 && *start) {
        fields[fcount++] = start;
    }
    if (fcount < 6) return false; // too few fields

    // Fields indices (1-based in spec): 1 QNAME,2 FLAG,3 RNAME,4 POS,5 MAPQ,6 CIGAR,7 RNEXT,8 PNEXT,9 TLEN,10 SEQ
    // Our array is 0-based, so:
    // flag -> fields[1], pos -> fields[3], mapq -> fields[4], cigar -> fields[5], seq -> fields[9]
    char *flag_str = (fcount > 1) ? fields[1] : NULL;
    char *pos_str = (fcount > 3) ? fields[3] : NULL;
    char *mapq_str = (fcount > 4) ? fields[4] : NULL;
    char *cigar = (fcount > 5) ? fields[5] : NULL;
    char *seq = (fcount > 9) ? fields[9] : NULL;

    if (!flag_str || !pos_str || !mapq_str || !cigar) return false;
    int flag = atoi(flag_str);
    int pos1 = atoi(pos_str); // SAM pos is 1-based; pos1==0 is possible for unmapped
    int mapq = atoi(mapq_str);
    int8_t strand = (flag & 0x10) ? -1 : 1;

    // if unmapped or rname "*" or cigar "*" skip
    if (pos1 <= 0 || cigar[0] == '*' || strcmp(cigar, "*") == 0) return false;

    uint32_t ref_pos0 = (uint32_t)(pos1 - 1);
    uint32_t query_pos0 = 0;
    uint32_t qlen = 0;
    if (seq && seq[0] != '*') {
        qlen = (uint32_t)strlen(seq);
    } else {
        // SEQ not available; we'll still parse cigar to get query advance but cannot validate length
        // Query length will be inferred from maximum query_end observed
    }

    // Before advancing cigar, set starting query_pos0 to 0; advance_cigar_blocks will handle leading S/I etc.
    size_t before_n = plot->n_segments;
    advance_cigar_blocks(cigar, &ref_pos0, &query_pos0, plot, strand, (uint8_t)mapq);

    // update maxima: ref_pos0 and query_pos0 are advanced to end of read (after processing cigar)
    if (ref_pos0 > *out_max_ref) *out_max_ref = ref_pos0;
    if (qlen > *out_max_query) {
        // if SEQ provided, consider full read length as potential max
        if (qlen > *out_max_query) *out_max_query = qlen;
    } else {
        if (query_pos0 > *out_max_query) *out_max_query = query_pos0;
    }

    // if append failed due to overflow, detect by comparing n_segments
    if (plot->n_segments == before_n) {
        // no segments appended (or overflow). If overflow, we stop further parsing
        return false;
    }
    return true;
}

// 主入口：从 SAM 文件构建 dotplot JSON，并写入 out_json_path（例如 "/dotplot.json"）
// 返回 true 表示成功
bool build_dotplot_json_from_sam(const char *sam_path, const char *out_json_path) {
    if (!sam_path || !out_json_path) return false;
    FILE *f = fopen(sam_path, "rb");
    if (!f) return false;
    dotplot_t plot;
    plot.n_segments = 0;
    plot.ref_length = 0;
    plot.query_length = 0;

    char *linebuf = (char*)malloc(SAM_LINE_BUF);
    if (!linebuf) { fclose(f); return false; }

    uint32_t max_ref = 0;
    uint32_t max_query = 0;

    while (fgets(linebuf, SAM_LINE_BUF, f)) {
        if (linebuf[0] == '@') {
            // header line - optionally parse @SQ for reference lengths (not required)
            // Example: @SQ SN:chr1 LN:248956422
            // We skip for now; using max_ref from alignments is sufficient
            continue;
        }
        // parse alignment line (we will modify linebuf in-place)
        // If parse fails due to segment overflow, we stop parsing but still write JSON of collected segments
        bool ok = parse_sam_alignment_line(linebuf, &plot, &max_ref, &max_query);
        if (!ok) {
            // if parse_sam_alignment_line returns false due to no segments (unmapped or invalid), continue
            // but if segments array is full, we should break
            if (plot.n_segments >= MAX_SEGMENTS) {
                // reached capacity; stop parsing more lines
                break;
            }
            // else continue scanning other lines
            continue;
        }
    }

    free(linebuf);
    fclose(f);

    // set lengths as maxima observed; ensure at least 1 to avoid divide-by-zero
    plot.ref_length = (max_ref > 0) ? max_ref : 1;
    plot.query_length = (max_query > 0) ? max_query : 1;

    // write to JSON
    bool wrote = dotplot_to_json_file(&plot, out_json_path);
    return wrote;
}


static int is_existing_file(const char *path) {
    if (!path) return 0;
    if (path[0] != '/') return 0;
    struct stat st;
    return (stat(path, &st) == 0);
}

// 最大允许直接传入内容的长度（若你仍支持通过字符串传内容）
#define MAX_INPUT_BYTES (100 * 1024 * 1024UL) // 100 MB，可按需调整

EMSCRIPTEN_KEEPALIVE
const char* run_minimap2_wasm(const char* args_str, const char* ref_content_or_path, const char* query_content_or_path) {
    if (!args_str || !ref_content_or_path || !query_content_or_path) {
        return strdup("Error: Null pointer input");
    }

    // 解析 args_str -> argv（两遍：先计数再构造）
    char *args_copy = strdup(args_str);
    if (!args_copy) return strdup("Error: Failed to duplicate args_str");
    int tmp_argc = 0;
    char *t = strtok(args_copy, " ");
    while (t) { tmp_argc++; t = strtok(NULL, " "); }
    free(args_copy);

    char **argv = (char**)malloc(sizeof(char*) * (tmp_argc + 4));
    if (!argv) return strdup("Error: Failed to allocate argv");

    args_copy = strdup(args_str);
    if (!args_copy) { free(argv); return strdup("Error: Failed to duplicate args_str"); }

    int argc = 0;
    argv[argc++] = "minimap2";
    t = strtok(args_copy, " ");
    while (t) {
        argv[argc++] = t;
        t = strtok(NULL, " ");
    }
    argv[argc] = NULL; // 必须以 NULL 结尾

    // 决定参考/查询文件路径：优先使用已存在的虚拟 FS 文件，否则把内容写入临时文件
    char *ref_filename = NULL;
    char *query_filename = NULL;
    int ref_is_temp = 0;
    int query_is_temp = 0;

    // helper 用于构建临时文件名
    long ts = (long)time(NULL);
    int pid = (int)getpid();
    char tmpbuf[128];

    // 处理 ref
    if (is_existing_file(ref_content_or_path)) {
        ref_filename = strdup(ref_content_or_path);
        if (!ref_filename) goto ERR_ALLOC;
    } else {
        size_t ref_len = strlen(ref_content_or_path);
        if (ref_len == 0 || ref_len > MAX_INPUT_BYTES) {
            goto ERR_INVALID_INPUT;
        }
        // 生成相对唯一临时路径
        snprintf(tmpbuf, sizeof(tmpbuf), "/target_%ld_%d.fa", ts, pid);
        FILE *fr = fopen(tmpbuf, "wb");
        if (!fr) goto ERR_FILE;
        if (fwrite(ref_content_or_path, 1, ref_len, fr) != ref_len) { fclose(fr); goto ERR_FILE; }
        fclose(fr);
        ref_filename = strdup(tmpbuf);
        if (!ref_filename) goto ERR_ALLOC;
        ref_is_temp = 1;
    }

    // 处理 query
    if (is_existing_file(query_content_or_path)) {
        query_filename = strdup(query_content_or_path);
        if (!query_filename) goto ERR_ALLOC;
    } else {
        size_t qlen = strlen(query_content_or_path);
        if (qlen == 0 || qlen > MAX_INPUT_BYTES) {
            goto ERR_INVALID_INPUT;
        }
        snprintf(tmpbuf, sizeof(tmpbuf), "/query_%ld_%d.fa", ts, pid);
        FILE *fq = fopen(tmpbuf, "wb");
        if (!fq) goto ERR_FILE;
        if (fwrite(query_content_or_path, 1, qlen, fq) != qlen) { fclose(fq); goto ERR_FILE; }
        fclose(fq);
        query_filename = strdup(tmpbuf);
        if (!query_filename) goto ERR_ALLOC;
        query_is_temp = 1;
    }

    // 生成唯一输出文件名（.sam 或 .log，根据需要）
    char stdout_path[128];
    char stderr_path[128];
    snprintf(stdout_path, sizeof(stdout_path), "/output_stdout_%ld_%d.sam", ts, pid);
    snprintf(stderr_path, sizeof(stderr_path), "/output_stderr_%ld_%d.log", ts, pid);
    char *output_stdout = strdup(stdout_path);
    if (!output_stdout) goto ERR_ALLOC;
    char *output_stderr = strdup(stderr_path);
    if (!output_stderr) { free(output_stdout); goto ERR_ALLOC; }

    // 保存原始 stdout/stderr
    int orig_stdout = dup(STDOUT_FILENO);
    int orig_stderr = dup(STDERR_FILENO);
    if (orig_stdout == -1 || orig_stderr == -1) {
        // 无法 dup，准备返回错误
        free(output_stdout);
        goto ERR_FD;
    }

    // 打开 stdout/stderr 文件并重定向
    FILE *fout = fopen(output_stdout, "wb+");
    if (!fout) {
        free(output_stdout); free(output_stderr);
        goto ERR_FD;
    }
    FILE *ferr = fopen(output_stderr, "wb+");
    if (!ferr) {
        fclose(fout);
        free(output_stdout); free(output_stderr);
        goto ERR_FD;
    }

    if (dup2(fileno(fout), STDOUT_FILENO) == -1) {
        fclose(fout); fclose(ferr);
        free(output_stdout); free(output_stderr);
        goto ERR_FD;
    }
    if (dup2(fileno(ferr), STDERR_FILENO) == -1) {
        // 恢复 stdout
        dup2(orig_stdout, STDOUT_FILENO);
        fclose(fout); fclose(ferr);
        free(output_stdout); free(output_stderr);
        goto ERR_FD;
    }
    // 关闭 FILE*，fd 已被 dup2 拷贝
    fclose(fout);
    fclose(ferr);

    // 运行 minimap2（注意：argv 应包含命令行参数，main.js 通常会把路径放到 args_str）
    int rc = minimap2_main(argc, argv);

    // 恢复 stdout/stderr
    fflush(stdout); fflush(stderr);
    dup2(orig_stdout, STDOUT_FILENO);
    dup2(orig_stderr, STDERR_FILENO);
    close(orig_stdout); close(orig_stderr);

    // 清理 argv 相关内存（args_copy 指向 argv 中的字符串）
    free(args_copy);
    free(argv);

    //
    // 生成 dotplot.json 供前端直接读取（方案 A）
    // 注意：出于效率考虑，你可以在这里忽略错误（仅记录），不影响返回 SAM 路径
    const char *dotpath = "/dotplot.json";
    if (!build_dotplot_json_from_sam(output_stdout, dotpath)) {
        // 可用 stderr 输出以便调试
        fprintf(stderr, "Warning: failed to build dotplot JSON from SAM %s\n", output_stdout);
    } else {
        // 成功写入 /dotplot.json
        fprintf(stderr, "WROTE %s\n", dotpath);
    }

    // 如果我们在 C 端创建了临时输入文件，删除它们以释放空间（输出文件保留给 JS）
    if (ref_is_temp && ref_filename) {
        unlink(ref_filename);
        free(ref_filename);
    } else if (ref_filename) {
        free(ref_filename);
    }

    if (query_is_temp && query_filename) {
        unlink(query_filename);
        free(query_filename);
    } else if (query_filename) {
        free(query_filename);
    }

    // 注意：我们返回 stdout 文件路径；前端若需要 stderr，可以按需读取 output_stderr 路径
    // 目前遵循你的要求：前端只在 UI 输出 stdout
    //free(output_stderr); // 若不想把 stderr 路径暴露给前端，可在此释放（若想保留，改为返回或传回）
    free(output_stdout);

    // 返回输出文件路径（由调用方通过 wasm_free 释放字符串，并读取/删除该虚拟文件）
    //return output_stdout;
    return output_stderr;

// ----- 错误处理 -----
ERR_FD:
    // 尝试恢复 stdout/stderr（若有）
    if (orig_stdout != -1) { dup2(orig_stdout, STDOUT_FILENO); close(orig_stdout); }
    if (orig_stderr != -1) { dup2(orig_stderr, STDERR_FILENO); close(orig_stderr); }
    if (args_copy) free(args_copy);
    if (argv) free(argv);
    if (ref_filename) { if (ref_is_temp) unlink(ref_filename); free(ref_filename); }
    if (query_filename) { if (query_is_temp) unlink(query_filename); free(query_filename); }
    return strdup("Error: Failed to setup file descriptors or output file");

ERR_FILE:
    if (args_copy) free(args_copy);
    if (argv) free(argv);
    if (ref_filename) { if (ref_is_temp) unlink(ref_filename); free(ref_filename); }
    if (query_filename) { if (query_is_temp) unlink(query_filename); free(query_filename); }
    return strdup("Error: Failed to write temporary input file");

ERR_ALLOC:
    if (args_copy) free(args_copy);
    if (argv) free(argv);
    if (ref_filename) { if (ref_is_temp) unlink(ref_filename); free(ref_filename); }
    if (query_filename) { if (query_is_temp) unlink(query_filename); free(query_filename); }
    return strdup("Error: Memory allocation failure");

ERR_INVALID_INPUT:
    if (args_copy) free(args_copy);
    if (argv) free(argv);
    if (ref_filename) { if (ref_is_temp) unlink(ref_filename); free(ref_filename); }
    if (query_filename) { if (query_is_temp) unlink(query_filename); free(query_filename); }
    return strdup("Error: Invalid or empty input content");
}
// 导出一个释放函数，JS 可调用 wasm_free(ptr) 来释放 C 端 strdup 的内存
EMSCRIPTEN_KEEPALIVE
void wasm_free(void *p) { if (p) free(p); }

int minimap2_main(int argc, char *argv[]) {
    // fprintf(stdout, "wasm_wrap.c(minimap2_main) 01 entered\n");//
    
    const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:b:O:E:m:N:Qu:R:hF:LC:yYPo:e:U:J:j:";
    ketopt_t o = KETOPT_INIT;
    mm_mapopt_t opt;
    mm_idxopt_t ipt;
    int i, c, n_threads = 1, n_parts, old_best_n = -1;
    float spsc_scale = 0.7f;
    char *fnw = 0, *rg = 0, *fn_bed_junc = 0, *fn_bed_jump = 0, *fn_bed_pass1 = 0, *fn_spsc = 0, *s, *alt_list = 0;
    FILE *fp_help = stderr;
    mm_idx_reader_t *idx_rdr;
    mm_idx_t *mi;
    
    // 关键修复：初始化 opt 结构体
    memset(&opt, 0, sizeof(mm_mapopt_t));
    // 强制禁用多线程
    opt.flag = 0;
    opt.flag &= ~MM_F_2_IO_THREADS;
    
    mm_verbose = 3;
    liftrlimit();
    mm_realtime0 = realtime();
    mm_set_opt(0, &ipt, &opt);
    
    // 第一轮参数解析
    while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
        // fprintf(stdout, "wasm_wrap.c(minimap2_main) 02 entered  option: '%c' (%d) ", c, c);//
        if (o.arg) {
            fprintf(stdout, "        ... with value: \"%s\"\n", o.arg);
        } else {
            // fprintf(stdout, "\n");//
        }
        
        if (c == 'x') {
            if (mm_set_opt(o.arg, &ipt, &opt) < 0) {
                fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
                return 1;
            }
        } else if (c == ':') {
            fprintf(stderr, "[ERROR] missing option argument\n");
            return 1;
        } else if (c == '?') {
            fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
            return 1;
        }
    }
    o = KETOPT_INIT;
    // fprintf(stdout, "wasm_wrap.c(minimap2_main) 03 before ketopt\n");//
    
    // 第二轮参数解析
    while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
        // fprintf(stdout, "wasm_wrap.c(minimap2_main) 04 ketopt found option: '%c' (ASCII: %d) ", c, c);//
        if (o.arg) {
            fprintf(stdout, "        ... with value: \"%s\"\n", o.arg);
        } else {
            // fprintf(stdout, "\n");//
        }
        
        if (c == 'w') ipt.w = atoi(o.arg);
        else if (c == 'k') ipt.k = atoi(o.arg);
        else if (c == 'H') ipt.flag |= MM_I_HPC;
        else if (c == 'd') fnw = o.arg;
        else if (c == 't') n_threads = atoi(o.arg);
        else if (c == 'v') mm_verbose = atoi(o.arg);
        else if (c == 'g') opt.max_gap = (int)mm_parse_num(o.arg);
        else if (c == 'G') mm_mapopt_max_intron_len(&opt, (int)mm_parse_num(o.arg));
        else if (c == 'F') opt.max_frag_len = (int)mm_parse_num(o.arg);
        else if (c == 'N') old_best_n = opt.best_n, opt.best_n = atoi(o.arg);
        else if (c == 'p') opt.pri_ratio = atof(o.arg);
        else if (c == 'M') opt.mask_level = atof(o.arg);
        else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
        else if (c == 'D') opt.flag |= MM_F_NO_DIAG;
        else if (c == 'P') opt.flag |= MM_F_ALL_CHAINS;
        else if (c == 'X') opt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
        else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
        else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
        else if (c == 'Y') opt.flag |= MM_F_SOFTCLIP;
        else if (c == 'L') opt.flag |= MM_F_LONG_CIGAR;
        else if (c == 'y') opt.flag |= MM_F_COPY_COMMENT;
        else if (c == 'T') opt.sdust_thres = atoi(o.arg);
        else if (c == 'n') opt.min_cnt = atoi(o.arg);
        else if (c == 'm') opt.min_chain_score = atoi(o.arg);
        else if (c == 'A') opt.a = atoi(o.arg);
        else if (c == 'B') opt.b = atoi(o.arg);
        else if (c == 'b') opt.transition = atoi(o.arg);
        else if (c == 's') opt.min_dp_max = atoi(o.arg);
        else if (c == 'C') opt.noncan = atoi(o.arg);
        else if (c == 'I') ipt.batch_size = mm_parse_num(o.arg);
        else if (c == 'K') opt.mini_batch_size = mm_parse_num(o.arg);
        else if (c == 'e') opt.occ_dist = mm_parse_num(o.arg);
        else if (c == 'R') rg = o.arg;
        else if (c == 'h') fp_help = stdout;
        else if (c == '2') opt.flag |= MM_F_2_IO_THREADS;
        else if (c == 'j') fn_bed_jump = o.arg;
        else if (c == 'J') {
            int t;
            t = atoi(o.arg);
            if (t == 0) opt.flag |= MM_F_SPLICE_OLD;
            else if (t == 1) opt.flag &= ~MM_F_SPLICE_OLD;
        } else if (c == 'o') {
            // fprintf(stdout, "wasm_wrap.c(minimap2_main) 05 check c==o o.arg==%s\n",o.arg);//
            // fprintf(stdout, "wasm_wrap.c(minimap2_main) 06 check freopen o.arg \n");//
            if (freopen(o.arg, "wb", stdout) == NULL) {
                fprintf(stderr, "[ERROR] failed to write the output to file  %s : %s\n", o.arg, strerror(errno));
                return 1;
            }
        } else if (c == 300) ipt.bucket_bits = atoi(o.arg);
        else if (c == 302) opt.seed = atoi(o.arg);
        else if (c == 303) mm_dbg_flag |= MM_DBG_NO_KALLOC;
        else if (c == 304) mm_dbg_flag |= MM_DBG_PRINT_QNAME;
        else if (c == 306) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, n_threads = 1;
        else if (c == 307) opt.max_chain_skip = atoi(o.arg);
        else if (c == 339) opt.max_chain_iter = atoi(o.arg);
        else if (c == 308) opt.min_ksw_len = atoi(o.arg);
        else if (c == 309) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, n_threads = 1;
        else if (c == 310) opt.flag |= MM_F_SPLICE;
        else if (c == 312) opt.flag |= MM_F_NO_LJOIN;
        else if (c == 317) opt.end_bonus = atoi(o.arg);
        else if (c == 318) opt.flag |= MM_F_INDEPEND_SEG;
        else if (c == 320) ipt.flag |= MM_I_NO_SEQ;
        else if (c == 321) opt.anchor_ext_shift = atoi(o.arg);
        else if (c == 322) opt.flag |= MM_F_FOR_ONLY;
        else if (c == 323) opt.flag |= MM_F_REV_ONLY;
        else if (c == 327) opt.max_clip_ratio = atof(o.arg);
        else if (c == 328) opt.min_mid_occ = atoi(o.arg);
        else if (c == 329) opt.flag |= MM_F_OUT_MD;
        else if (c == 331) opt.sc_ambi = atoi(o.arg);
        else if (c == 332) opt.flag |= MM_F_EQX;
        else if (c == 333) opt.flag |= MM_F_PAF_NO_HIT;
        else if (c == 334) opt.split_prefix = o.arg;
        else if (c == 335) opt.flag |= MM_F_NO_END_FLT;
        else if (c == 336) opt.flag |= MM_F_HARD_MLEVEL;
        else if (c == 337) opt.max_sw_mat = mm_parse_num(o.arg);
        else if (c == 338) opt.max_qlen = mm_parse_num(o.arg);
        else if (c == 340) fn_bed_junc = o.arg;
        else if (c == 341) opt.junc_bonus = atoi(o.arg);
        else if (c == 342) opt.flag |= MM_F_SAM_HIT_ONLY;
        else if (c == 343) opt.chain_gap_scale = atof(o.arg);
        else if (c == 351) opt.chain_skip_scale = atof(o.arg);
        else if (c == 344) alt_list = o.arg;
        else if (c == 345) opt.alt_drop = atof(o.arg);
        else if (c == 346) opt.mask_len = mm_parse_num(o.arg);
        else if (c == 348) opt.flag |= MM_F_QSTRAND | MM_F_NO_INV;
        else if (c == 349) opt.cap_kalloc = mm_parse_num(o.arg);
        else if (c == 350) opt.q_occ_frac = atof(o.arg);
        else if (c == 352) mm_dbg_flag |= MM_DBG_PRINT_CHAIN;
        else if (c == 353) opt.flag |= MM_F_NO_HASH_NAME;
        else if (c == 354) opt.flag |= MM_F_SECONDARY_SEQ;
        else if (c == 355) opt.flag |= MM_F_OUT_DS;
        else if (c == 356) opt.rmq_inner_dist = mm_parse_num(o.arg);
        else if (c == 357) fn_spsc = o.arg;
        else if (c == 360) opt.jump_min_match = mm_parse_num(o.arg);
        else if (c == 361) opt.flag |= MM_F_OUT_JUNC | MM_F_CIGAR;
        else if (c == 362) fn_bed_pass1 = o.arg;
        else if (c == 501) mm_dbg_flag |= MM_DBG_SEED_FREQ;
        else if (c == 363) spsc_scale = atof(o.arg);
        else if (c == 358 || c == 364) opt.junc_pen = atoi(o.arg);
        else if (c == 330) {
            fprintf(stderr, "[WARNING] \033[1;31m --lj-min-ratio has been deprecated.\033[0m\n");
        } else if (c == 313) {
            if (o.arg == 0 || strcmp(o.arg, "dna") == 0) {
                opt.flag |= MM_F_SR;
            } else if (strcmp(o.arg, "rna") == 0) {
                opt.flag |= MM_F_SR_RNA;
            } else if (strcmp(o.arg, "no") == 0) {
                opt.flag &= ~(uint64_t)(MM_F_SR|MM_F_SR_RNA);
            } else if (mm_verbose >= 2) {
                opt.flag |= MM_F_SR;
                fprintf(stderr, "[WARNING]\033[1;31m --sr only takes 'dna' or 'rna'. Invalid values are assumed to be 'dna'.\033[0m\n");
            }
        } else if (c == 314) {
            yes_or_no(&opt, MM_F_FRAG_MODE, o.longidx, o.arg, 1);
        } else if (c == 315) {
            yes_or_no(&opt, MM_F_NO_PRINT_2ND, o.longidx, o.arg, 0);
        } else if (c == 316) {
            opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
            if (o.arg == 0 || strcmp(o.arg, "short") == 0) {
                opt.flag &= ~MM_F_OUT_CS_LONG;
            } else if (strcmp(o.arg, "long") == 0) {
                opt.flag |= MM_F_OUT_CS_LONG;
            } else if (strcmp(o.arg, "none") == 0) {
                opt.flag &= ~MM_F_OUT_CS;
            } else if (mm_verbose >= 2) {
                fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
            }
        } else if (c == 319) {
            yes_or_no(&opt, MM_F_SPLICE_FLANK, o.longidx, o.arg, 1);
        } else if (c == 324) {
            yes_or_no(&opt, MM_F_HEAP_SORT, o.longidx, o.arg, 1);
        } else if (c == 326) {
            yes_or_no(&opt, MM_F_NO_DUAL, o.longidx, o.arg, 0);
        } else if (c == 347) {
            if (o.arg) yes_or_no(&opt, MM_F_RMQ, o.longidx, o.arg, 1);
            else opt.flag |= MM_F_RMQ;
        } else if (c == 359) {
            if (strcmp(o.arg, "no") == 0) opt.flag |= MM_F_INDEPEND_SEG;
            else if (strcmp(o.arg, "weak") == 0) opt.flag |= MM_F_WEAK_PAIRING, opt.flag &= ~(uint64_t)MM_F_INDEPEND_SEG;
            else {
                if (strcmp(o.arg, "strong") != 0 && mm_verbose >= 2)
                    fprintf(stderr, "[WARNING]\033[1;31m unrecognized argument for --pairing; assuming 'strong'.\033[0m\n");
                opt.flag &= ~(uint64_t)(MM_F_INDEPEND_SEG|MM_F_WEAK_PAIRING);
            }
        } else if (c == 'S') {
            opt.flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
            if (mm_verbose >= 2)
                fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
        } else if (c == 'V') {
            puts(MM_VERSION);
            return 0;
        } else if (c == 'r') {
            opt.bw = (int)mm_parse_num2(o.arg, &s);
            if (*s == ',') opt.bw_long = (int)mm_parse_num2(s + 1, &s);
        } else if (c == 'U') {
            opt.min_mid_occ = strtol(o.arg, &s, 10);
            if (*s == ',') opt.max_mid_occ = strtol(s + 1, &s, 10);
        } else if (c == 'f') {
            double x;
            char *p;
            x = strtod(o.arg, &p);
            if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
            else opt.mid_occ = (int)(x + .499);
            if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
        } else if (c == 'u') {
            if (*o.arg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV;
            else if (*o.arg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV;
            else if (*o.arg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR;
            else if (*o.arg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV);
            else {
                fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
                return 1;
            }
        } else if (c == 'z') {
            opt.zdrop = opt.zdrop_inv = strtol(o.arg, &s, 10);
            if (*s == ',') opt.zdrop_inv = strtol(s + 1, &s, 10);
        } else if (c == 'O') {
            opt.q = opt.q2 = strtol(o.arg, &s, 10);
            if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
        } else if (c == 'E') {
            opt.e = opt.e2 = strtol(o.arg, &s, 10);
            if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
        }
    }
    fprintf(stderr, "wasm_wrap.c(minimap2_main) 02 set options\n");
    
    // 参数校验
    if (!fnw && !(opt.flag&MM_F_CIGAR))
        ipt.flag |= MM_I_NO_SEQ;
    
    if (mm_check_opt(&ipt, &opt) < 0)
        return 1;
    
    if (opt.best_n == 0) {
        fprintf(stderr, "[WARNING]\033[1;31m changed '-N 0' to '-N %d --secondary=no'.\033[0m\n", old_best_n);
        opt.best_n = old_best_n, opt.flag |= MM_F_NO_PRINT_2ND;
    }
    
    // 帮助信息
    if (argc == o.ind || fp_help == stdout) {
        fprintf(fp_help, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "  Indexing:\n");
        fprintf(fp_help, "    -H           use homopolymer-compressed k-mer (preferrable for PacBio)\n");
        fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", ipt.k);
        fprintf(fp_help, "    -w INT       minimizer window size [%d]\n", ipt.w);
        fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [8G]\n");
        fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
        fprintf(fp_help, "  Mapping:\n");
        fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
        fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
        fprintf(fp_help, "    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
        fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
        fprintf(fp_help, "    -r NUM[,NUM] chaining/alignment bandwidth and long-join bandwidth [%d,%d]\n", opt.bw, opt.bw_long);
        fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
        fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
        fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
        fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
        fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
        fprintf(fp_help, "  Alignment:\n");
        fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
        fprintf(fp_help, "    -B INT       mismatch penalty (larger value for lower divergence) [%d]\n", opt.b);
        fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
        fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
        fprintf(fp_help, "    -z INT[,INT] Z-drop score and inversion Z-drop score [%d,%d]\n", opt.zdrop, opt.zdrop_inv);
        fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
        fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
        fprintf(fp_help, "    -J INT       splice mode. 0: original minimap2 model; 1: miniprot model [1]\n");
        fprintf(fp_help, "    -j FILE      junctions in BED12 to extend *short* RNA-seq alignment []\n");
        fprintf(fp_help, "  Input/Output:\n");
        fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
        fprintf(fp_help, "    -o FILE      output alignments to FILE [stdout]\n");
        fprintf(fp_help, "    -L           write CIGAR with >65535 ops at the CG tag\n");
        fprintf(fp_help, "    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
        fprintf(fp_help, "    -c           output CIGAR in PAF\n");
        fprintf(fp_help, "    --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]\n");
        fprintf(fp_help, "    --ds         output the ds tag, which is an extension to cs\n");
        fprintf(fp_help, "    --MD         output the MD tag\n");
        fprintf(fp_help, "    --eqx        write =/X CIGAR operators\n");
        fprintf(fp_help, "    -Y           use soft clipping for supplementary alignments\n");
        fprintf(fp_help, "    -y           copy FASTA/Q comments to output SAM\n");
        fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
        fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
        fprintf(fp_help, "    --version    show version number\n");
        fprintf(fp_help, "  Preset:\n");
        fprintf(fp_help, "    -x STR       preset (always applied before other options; see minimap2.1 for details) []\n");
        fprintf(fp_help, "                 - lr:hq - accurate long reads (error rate <1%%) against a reference genome\n");
        fprintf(fp_help, "                 - splice/splice:hq - spliced alignment for long reads/accurate long reads\n");
        fprintf(fp_help, "                 - splice:sr - spliced alignment for short RNA-seq reads\n");
        fprintf(fp_help, "                 - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5%% sequence divergence\n");
        fprintf(fp_help, "                 - sr - short reads against a reference\n");
        fprintf(fp_help, "                 - map-pb/map-hifi/map-ont/map-iclr - CLR/HiFi/Nanopore/ICLR vs reference mapping\n");
        fprintf(fp_help, "                 - ava-pb/ava-ont - PacBio CLR/Nanopore read overlap\n");
        fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of these and other advanced command-line options.\n");
        return fp_help == stdout? 0 : 1;
    }
    fprintf(stderr, "wasm_wrap.c(minimap2_main) 03 set options\n");
    print_minimap2_options(&ipt, &opt);
    
    // 输入校验
    if ((opt.flag & MM_F_SR) && argc - o.ind > 3) {
        fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
        return 1;
    }
    
    fprintf(stderr, "wasm_wrap.c(minimap2_main) 04 read argv[o.ind] ref/target file %s \n", argv[o.ind]);
    
    // 打开索引/参考文件
    idx_rdr = mm_idx_reader_open(argv[o.ind], &ipt, fnw);
    if (idx_rdr == 0) {
        fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
        return 1;
    }
    fprintf(stderr, "wasm_wrap.c(minimap2_main) 05 opened ref/target file\n");
    
    // 校验
    if (!idx_rdr->is_idx && fnw == 0 && argc - o.ind < 2) {
        fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
        mm_idx_reader_close(idx_rdr);
        return 1;
    }
    
    fprintf(stderr, "wasm_wrap.c(minimap2_main) 06 mm_verbose\n");
    
    if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
        fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");
    
    fprintf(stderr, "wasm_wrap.c(minimap2_main) 07 mm_idx_reader_read idx_rdr\n");
    
    // 核心比对逻辑
    while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
        fprintf(stderr, "wasm_wrap.c(minimap2_main) %s read %u sequences\n", argv[o.ind], mi->n_seq);
        
        if (mi->n_seq == 0) {
            mm_idx_destroy(mi);
            break;
        }
        
        long long total_bases_in_chunk = 0;
        for (int i = 0; i < mi->n_seq; ++i) {
            total_bases_in_chunk += mi->seq[i].len;
        }
        fprintf(stderr, "[INFO] 本次读取的数据块中包含 %d 条序列, 共计 %lld 个碱基。\n", mi->n_seq, total_bases_in_chunk);
        
        int ret;
        
        if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
            fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
            mm_idx_destroy(mi);
            mm_idx_reader_close(idx_rdr);
            return 1;
        }
        
        if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
            if (mm_idx_reader_eof(idx_rdr)) {
                if (opt.split_prefix == 0)
                    ret = mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
                else
                    ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
            } else {
                ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
                if (opt.split_prefix == 0 && mm_verbose >= 2)
                    fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.\033[0m\n");
            }
            if (ret != 0) {
                mm_idx_destroy(mi);
                mm_idx_reader_close(idx_rdr);
                return 1;
            }
        }
        
        if (mm_verbose >= 3)
            fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
                    __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
        
        if (argc != o.ind + 1) mm_mapopt_update(&opt, mi);
        
        if (mm_verbose >= 3) mm_idx_stat(mi);
        
        // 加载剪接位点
        if (fn_bed_junc) {
            mm_idx_bed_read(mi, fn_bed_junc, 1);
            if (mi->I == 0 && mm_verbose >= 2)
                fprintf(stderr, "[WARNING] failed to load the junction BED file\n");
        }
        if (fn_bed_jump) {
            mm_idx_jjump_read(mi, fn_bed_jump, MM_JUNC_ANNO, -1);
            if (mi->J == 0 && mm_verbose >= 2)
                fprintf(stderr, "[WARNING] failed to load the jump BED file\n");
        }
        if (fn_bed_pass1) {
            mm_idx_jjump_read(mi, fn_bed_pass1, MM_JUNC_MISC, 5);
            if (mi->J == 0 && mm_verbose >= 2)
                fprintf(stderr, "[WARNING] failed to load the pass-1 jump BED file\n");
        }
        if (fn_spsc) {
            mm_idx_spsc_read2(mi, fn_spsc, mm_max_spsc_bonus(&opt), spsc_scale);
            if (mi->spsc == 0 && mm_verbose >= 2)
                fprintf(stderr, "[WARNING] failed to load the splice score file\n");
        }
        if (alt_list) mm_idx_alt_read(mi, alt_list);
        
        // 比对逻辑
        if (argc - (o.ind + 1) == 0) {
            mm_idx_destroy(mi);
            continue;
        }
        ret = 0;
        
        if (!(opt.flag & MM_F_FRAG_MODE)) {
            for (i = o.ind + 1; i < argc; ++i) {
                fprintf(stderr, "read query file %s\n",argv[i]);
                ret = mm_map_file(mi, argv[i], &opt, n_threads);
                if (ret < 0) break;
            }
        } else {
            ret = mm_map_file_frag(mi, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_threads);
        }
        
        mm_idx_destroy(mi);
        if (ret < 0) {
            fprintf(stderr, "ERROR: failed to map the query file\n");
            return 1;
        }
    }
    
    n_parts = idx_rdr->n_parts;
    mm_idx_reader_close(idx_rdr);
    
    if (opt.split_prefix)
        mm_split_merge(argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_parts);
    
    if (fflush(stdout) == EOF) {
        perror("[ERROR] failed to write the results");
        return 1;
    }
    
    if (mm_verbose >= 3) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mm_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }
    
    return 0;
}