# IRSSG Unified Build

这是一个重新组织的IRSSG项目，将所有源码合并到一个统一的目录结构中，支持一次性编译。

## 目录结构

```
unified_irssg/
├── lib/           # 库源码文件 (原lib_ssg目录)
├── src/           # 主程序源码文件 (原src_irssg目录)
├── Makefile       # 统一编译文件
├── libIRSSG.a     # 生成的静态库
├── irssg          # 生成的可执行文件
└── README.md      # 本说明文件
```

## 编译方法

1. 加载Intel OneAPI环境：
   ```bash
   module load oneapi22.3
   ```

2. 进入项目目录：
   ```bash
   cd unified_irssg
   ```

3. 编译项目：
   ```bash
   make
   ```

4. 清理编译文件：
   ```bash
   make clean
   ```

## 编译流程

新的Makefile会自动：
1. 编译lib目录下的所有库源码文件
2. 创建静态库libIRSSG.a
3. 编译src目录下的所有主程序源码文件
4. 链接生成最终的可执行文件irssg

## 优势

- 简化了编译流程，从原来的两步编译（lib_ssg -> src_irssg）简化为一步
- 统一的目录结构，便于管理和维护
- 保持了原有的编译选项和依赖关系

## 单一可执行 + 运行时选择

已将两个主程序合并为一个可执行文件 `irssg`，通过参数在运行时选择模式：

- 默认执行普通（PW）路径
- 加 `--wann` 或 `--mode wann` 执行 Wannier 路径

示例：

```
./irssg -nk 1 10                      # 普通路径
./irssg --wann -nk 1 10               # Wannier 路径
./irssg --mode wann -tolE 1e-3        # Wannier 路径（等效）
```

说明：

- 其余原有参数不变；新增的 `--wann`/`--mode` 对旧代码是“未知参数”，内部会自动忽略，不影响参数解析。
- 源码层面对 `src_wann` 的 `comms/init/wave_data` 模块做了重命名（增加 `_wann` 后缀），以避免与普通路径的同名模块冲突，并将两套逻辑作为 `driver` 子程序由统一 `main` 调度。

## 集成 MOM2SSG（-ssg）

已将 Python 包 MOM2SSG 集成到安装包中，使用方式：

- 运行 Fortran 主流程（默认）:
  - `irssg [原有参数]`
- 调用 MOM2SSG（Python）:
  - `irssg -ssg [MOM2SSG参数]`

说明：
- `-ssg` 会转发后续参数给 `MOM2SSG.MOM2SSG:main`，等价于执行 `python -m MOM2SSG.MOM2SSG ...`。
- MOM2SSG 依赖 `numpy/scipy/spglib/pymatgen/phonopy`，请在运行 `-ssg` 前确保环境已安装这些依赖。

## pip 安装（本地构建）

支持通过 pip 本地安装，自动编译 Fortran 并安装命令 `irssg`：

环境要求：
- Intel ifort/ifx 编译器与 MKL（Makefile 使用 ifort + `-qmkl`）
- Python >= 3.8

本地安装（离线/集群环境建议禁用隔离以避免网络下载）：

```
python -m pip install . --no-build-isolation
```

说明：
- 构建阶段会执行 `make USE_ABS_DATA_PATH=0`，二进制与 `kLittleGroups` 数据将被打包进 Python 包中。
- 运行时 `irssg` 入口会自动设置 `IRVSPDATA` 指向包内数据目录，无需手动配置。
- 如需重新打包 wheel：`python -m pip wheel . --no-build-isolation`，生成的平台 wheel 将包含已编译的二进制。





# irssg_py2
