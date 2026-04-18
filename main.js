// 获取HTML元素
const outputElement = document.getElementById('output');
const refFileInput = document.getElementById('ref-file-input');
const queryFileInput = document.getElementById('query-file-input');
const runButton = document.getElementById('run-button');
const paramsInput = document.getElementById('params-input');

const refInfo = document.getElementById('ref-info');
const queryInfo = document.getElementById('query-info');

const downloadRstBtn = document.getElementById('download-button');

const canvas = document.getElementById('dotplot-canvas');
const ctx = canvas.getContext('2d');
const scoreFilter = document.getElementById('score-filter');
const scoreValue = document.getElementById('score-value');
const segmentList = document.getElementById('segment-list');

let resultPtr = null; // 初始化 resultPtr 为 null
let latestDownload = null; // { blob, name, url }

let latestDotplotData = null;

function renderDotplot(data, minScore = 0) {
    if (!data) return;
    latestDotplotData = data;

    const { ref_length, query_length, segments } = data;
    const width = canvas.width;
    const height = canvas.height;
    const xScale = width / Math.max(1, ref_length);
    const yScale = height / Math.max(1, query_length);

    ctx.clearRect(0, 0, width, height);

    // 画主对角线
    ctx.strokeStyle = '#bbb';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(width, height);
    ctx.stroke();

    const filtered = segments
        .filter(seg => seg.score >= minScore)
        .sort((a, b) => (b.ref_end - b.ref_start) - (a.ref_end - a.ref_start));

    filtered.forEach(seg => {
        const segWidth = Math.max(1, (seg.ref_end - seg.ref_start) * xScale);
        const segHeight = Math.max(1, (seg.query_end - seg.query_start) * yScale);
        const x = seg.ref_start * xScale;
        let y = seg.query_start * yScale;

        if (seg.strand < 0) {
            y = height - (seg.query_end * yScale);
        }

        ctx.fillStyle = seg.strand < 0
            ? 'rgba(220, 80, 80, 0.6)'
            : 'rgba(60, 140, 230, 0.6)';
        ctx.fillRect(x, y, segWidth, segHeight);
    });

    // 高分 segment 列表
    segmentList.innerHTML = filtered.slice(0, 5).map(seg => {
        return `<li>ref:[${seg.ref_start}-${seg.ref_end}], query:[${seg.query_start}-${seg.query_end}], score:${seg.score}, strand:${seg.strand}</li>`;
    }).join('');
}

// 假设你从 wasm 获取 JSON，调用下面函数
function onDotplotJsonReceived(jsonText) {
    try {
        const data = JSON.parse(jsonText);
        renderDotplot(data, Number(scoreFilter.value));
    } catch (err) {
        console.error("解析 dotplot JSON 失败", err);
    }
}

// helper: 把 File 写入 Emscripten FS（path 必须以 '/' 开头）
async function writeFileToFS(Module, file, path) {
    const ab = await file.arrayBuffer();
    const u8 = new Uint8Array(ab);
    try { Module.FS_unlink(path); } catch (e) { /* ignore */ }
    Module.FS_createDataFile('/', path.slice(1), u8, true, true);
}

// helper: 安全 unlink（忽略不存在的情况）
function safeUnlink(Module, path) {
    try { Module.FS_unlink(path); } catch (e) { /* ignore */ }
}

// 解码 Uint8Array -> UTF-8 字符串（适合 < 1MB）
function uint8ToString(u8) {
    return new TextDecoder().decode(u8);
}

// 下载
function cleanupDownload() {
    if (latestDownload && latestDownload.url) {
        URL.revokeObjectURL(latestDownload.url);
    }
    latestDownload = null;
    downloadRstBtn.disabled = true;
}
function prepareDownload(fileData, resultPath) {
    //const downloadName = resultPath.split('/').pop() || 'alignment.sam';
    const downloadName = 'output.sam';
    cleanupDownload(); // 先释放旧的
    const blob = new Blob([fileData], { type: 'application/octet-stream' });
    const url = URL.createObjectURL(blob);
    latestDownload = { blob, name: downloadName, url };
    downloadRstBtn.disabled = false;
}

// createMinimap2Module 是由你的 EXPORT_NAME 参数定义的全局函数。
createMinimap2Module().then(Module => {
    console.log("main.js: Minimap2 Wasm module loaded.");
    runButton.textContent = 'Run Alignment';
    runButton.disabled = false; 

    // ---- 使用 cwrap 创建 JS 友好型函数 ----
    const run_minimap2 = Module.cwrap(
        'run_minimap2_wasm', 
        'number',            
        ['string', 'string', 'string'] 
    );

    const free_result = Module.cwrap(
        'wasm_free',
        null,
        ['number']
    );

    // ---- 单击事件处理：写入 FS -> 调用 wasm -> 读取输出 -> 下载 -> 清理 ----
    async function onRunButtonClick() {
        cleanupDownload(); // 释放上次的结果
        const refFile = refFileInput.files[0];
        const queryFile = queryFileInput.files[0];

        if (!refFile || !queryFile) {
            outputElement.textContent = "Error: Please select both a reference and a query file.";
            return;
        }
        
        // 每次运行前重置并禁用按钮，防止并发
        resultPtr = null;
        runButton.disabled = true;
        outputElement.textContent = '';

        // 生成唯一的临时路径，避免冲突
        const ts = Date.now();
        const refPath = `/target_${ts}.fa`;
        const queryPath = `/query_${ts}.fa`;

        try {
            outputElement.textContent = "Writing files to WASM FS...";
            console.log(`main.js: Writing ${refFile.name} -> ${refPath} (${refFile.size} bytes)`);
            console.log(`main.js: Writing ${queryFile.name} -> ${queryPath} (${queryFile.size} bytes)`);

            await writeFileToFS(Module, refFile, refPath);
            await writeFileToFS(Module, queryFile, queryPath);

            outputElement.textContent = "Files written. Running alignment (this may take some time)...";
            const userParams = paramsInput.value.trim();
            const commandLineArgs = `${userParams} ${refPath} ${queryPath}`;
            console.log("main.js: Command:", commandLineArgs);

            // ---- 调用 Wasm 函数 ----
            resultPtr = Number(run_minimap2(commandLineArgs, refPath, queryPath));
            if (resultPtr === 0 || isNaN(resultPtr)) {
                outputElement.textContent = "Error: Wasm function returned a null pointer. The input might be invalid or an internal error occurred.";
                return;
            }

            // ---- 从 Wasm 内存中读取返回的文件路径（短字符串） ----
            const resultPath = Module.UTF8ToString(resultPtr);
            console.log("main.js: Result path:", resultPath);

            // ---- 以二进制方式读取输出文件（Uint8Array） ----
            let fileData;
            try {
                fileData = Module.FS_readFile(resultPath); // Uint8Array
            } catch (fsErr) {
                console.error("main.js: Failed to read output file from FS:", fsErr);
                outputElement.textContent = "Error: Failed to read alignment output file from virtual FS.";
                return;
            }
            const fileSize = fileData.length;
            console.log("main.js: Output file size:", fileSize, "bytes");

            // ---- 仅预览前 N 字节（避免展示巨大内容） ----
            const PREVIEW_BYTES = 2048; // 2 KB
            const previewSize = Math.min(PREVIEW_BYTES, fileData.length);
            const previewText = new TextDecoder().decode(fileData.subarray(0, previewSize));
            outputElement.textContent = previewText + `...\n\n Full file size: ${fileSize} bytes`;

            // ---- 触发下载（不自动）：准备 blob，但仅启用按钮 ----
            prepareDownload(fileData, resultPath);

            //
            // 读取 dotplot.json（方案 A：从虚拟 FS 读取）
            let dotJsonText = null;
            try {
                    // 优先尝试以 text 直接读取（若 Emscripten 版本支持 encoding 选项）
                    try {
                        dotJsonText = Module.FS_readFile('/dotplot.json', { encoding: 'utf8' });
                    } catch (e) {
                    // 如果上面不支持 encoding，退回到读取 Uint8Array 再解码
                    const dotArr = Module.FS_readFile('/dotplot.json'); // Uint8Array
                    if (dotArr && dotArr.length > 0) {
                        dotJsonText = new TextDecoder().decode(dotArr);
                    }
                    }
                    if (dotJsonText) {
                        onDotplotJsonReceived(dotJsonText);
                    } else {
                        console.warn('main.js: /dotplot.json not found or empty.');
                    }
                } catch (readErr) {
                    console.warn('main.js: Error reading /dotplot.json from WASM FS', readErr);
            }


            // ---- 清理：释放 C 端返回的 strdup 指针，再删除虚拟文件 ----
            try { free_result(resultPtr); } catch (freeErr) { console.warn("main.js: Failed to free C-side memory:", freeErr); }
            try { Module.FS_unlink(resultPath); } catch (unlinkErr) { console.warn("main.js: Failed to unlink virtual output file:", unlinkErr); }

        } catch (error) {
            console.error("main.js:  An error occurred:", error);
            outputElement.textContent = `An error occurred: ${error && error.message ? error.message : error}`;
        } finally {
            // 无论成功或失败，都删除临时输入文件，恢复按钮状态
            try { safeUnlink(Module, refPath); } catch(e) {}
            try { safeUnlink(Module, queryPath); } catch(e) {}
            runButton.disabled = false;
            resultPtr = null; 
        }
    }

    // 读取序列信息和长度
    async function readFastaInfo(file) {
        if (!file) return null;
        const text = await file.text();
        let name = '';
        let length = 0;
        for (const line of text.split(/\r?\n/)) {
            if (!line) continue;
            if (line.startsWith('>')) {
                if (!name) name = line.slice(1).trim() || file.name;
                continue;
            }
            length += line.trim().length;
        }
        return { name: name || file.name, length };
    }
    async function updateFileInfo(input, infoElement) {
        const info = await readFastaInfo(input.files[0]);
        infoElement.style.whiteSpace = 'pre-line';
        if (!info) { infoElement.textContent = '\n\n'; return; }
        infoElement.textContent = `${info.name} (${info.length} bp)\n`;
    }

    // 绑定事件
    refFileInput.addEventListener('change', () => updateFileInfo(refFileInput, refInfo));
    queryFileInput.addEventListener('change', () => updateFileInfo(queryFileInput, queryInfo));
    runButton.addEventListener('click', onRunButtonClick);
    downloadRstBtn.addEventListener('click', () => {
        if (!latestDownload) return;
        const a = document.createElement('a');
        a.href = latestDownload.url;
        a.download = latestDownload.name;
        document.body.appendChild(a);
        a.click();
        a.remove();
    });
    
    scoreFilter.addEventListener('input', () => {
        const minScore = Number(scoreFilter.value);
        scoreValue.textContent = minScore;
        renderDotplot(latestDotplotData, minScore);
    });

}).catch(err => {
    console.error("main.js: Failed to initialize wasm module:", err);
    if (outputElement) outputElement.textContent = "Failed to initialize wasm module. Check console for details.";
});