# Hypergraph Partitioner (hyperml)

## 概述

本項目實現了一個基於 hMETIS Heavy-Edge Matching (HEM) 的 2-way multilevel hypergraph partitioner。該算法包含以下四個主要階段：

1. **Coarsening**：使用 HEM 將原始超圖逐層粗化。
2. **Initial Partitioning**：在最粗化層級上，以 Greedy + FM 演算法產生初始分割。
3. **Projection**：將分割結果投影回下一層粗化前的超圖。
4. **Refinement**：在每層投影後使用超圖版 FM 精煉，並滿足 balance constraint 與最小化 cut size。

## 功能與需求

- **2-way partitioning** (k=2)。
- 結果需同時滿足：
  - **Balance constraint**：各 partition 的頂點總權重在允許偏差 ε 範圍內。
  - **最小化** hyperedge cut size。
- **Coarsening**：
  - 採用 Heavy-Edge Matching (HEM)。
  - 原始頂點 & 超邊權重均為單位權重 (1)；粗化後權重為合併後之總和。
  - 控制粗化聚合大小，避免過度不均衡。
- **Initial Partitioning**：
  - 處理非單位頂點/超邊權重。
  - 演算法：Greedy-fill 按頂點權重 + FM 調整。
- **Refinement**：
  - 採用 hypergraph FM 演算法，考慮多頂點移動增益與 balance。
- 支援超圖規模至百萬頂點以下。

## 使用方法

### 安裝依賴
本項目使用 C 語言實現，需要 GCC 編譯器。請確保已安裝 GCC。

### 編譯
在項目根目錄下執行以下命令進行編譯：
```sh
make
```

### 執行
編譯成功後，可使用以下命令運行程序：
```sh
./hyperml <input.hgr> [epsilon]
```
- `<input.hgr>`：輸入的超圖文件，格式為 `#nets #vertices`，後續每行列出 hyperedge。
- `[epsilon]`：可選參數，平衡容差值，默認為 0.03。

### 輸出
程序將輸出以下信息：
```
CutSize <cut_size>
Partition Sizes: <partition A size, partition B size>
Balance Deviation: <maximum balance deviation>
Total Execution Time: <total execution time>
```

## 資料結構

1. **Vertex**：
   - `struct V { int id; size_t weight; int match; }`；
   - `match=-1` 表示未匹配。

2. **Hyperedge**：
   - `struct E { vector<int> pins; size_t weight; }`.

3. **Hypergraph**：
   - `vector<V> vertices; vector<E> edges;`  
   - 建立「頂點→超邊」反向索引以加速匹配查詢。

4. **Coarsening Maps**：
   - `vector<int> coarse_map`：原頂點到粗化頂點 ID 的映射。
   - 新增粗化層級的 `Hypergraph` 實例。

## 參考資料

- hMETIS: [https://people.eecs.berkeley.edu/~smarkar/hmetis/](https://people.eecs.berkeley.edu/~smarkar/hmetis/)
```