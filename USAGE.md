# POS2SSG 使用说明

## 安装

### 方法1：使用安装脚本（推荐）
```bash
./install.sh
```

### 方法2：手动安装
```bash
# 构建包
python -m build

# 安装包
pip install dist/pos2ssg-1.0.0-py3-none-any.whl
```

## 使用方法

### 基本命令
```bash
# 查看帮助
POS2SSG --help

# 分析POSCAR文件，查找可能的SSG编号
POS2SSG -c POSCAR

# 查找特定Ik值的SSG编号
POS2SSG -c POSCAR -Ik 2

# 生成特定SSG的磁结构
POS2SSG -c POSCAR -ssg 43.2.4.1.P

# 使用自定义磁原子和磁矩
POS2SSG -c POSCAR -ssg 43.2.4.1.P --magatom Fe Cu --magmom 3.0 2.0

# 调整对称性分析精度
POS2SSG -c POSCAR -ssg 43.2.4.1.P --tolerance 1e-3

# 指定输出文件名
POS2SSG -c POSCAR -ssg 43.2.4.1.P --output my_structure.vasp

# 启用详细输出
POS2SSG -c POSCAR -ssg 43.2.4.1.P --verbose
```

### Python API使用
```python
from POS2SSG import main, construct_magmom_dict, read_poscar, write_poscar

# 读取结构
cell = read_poscar('POSCAR')

# 设置磁矩
magmom_dict = construct_magmom_dict(['Fe', 'Cu'], [3.0, 2.0])

# 使用其他功能...
```

## 文件结构

```
pos2ssg/
├── POS2SSG/              # 源代码目录
│   ├── __init__.py
│   ├── POS2SSG.py        # 主程序
│   ├── deal_cell.py      # 晶格处理
│   ├── generate_mag.py   # 磁矩生成
│   ├── generate_structure.py
│   ├── mag_atom.py       # 磁原子处理
│   └── poscar_io.py      # POSCAR文件IO
├── pyproject.toml        # 包配置
├── MANIFEST.in          # 打包清单
├── README.md            # 项目说明
├── LICENSE              # 许可证
├── requirements.txt     # 依赖列表
├── install.sh          # 安装脚本
└── dist/               # 构建产物
    ├── pos2ssg-1.0.0-py3-none-any.whl
    └── pos2ssg-1.0.0.tar.gz
```

## 依赖项

- Python >= 3.8
- numpy >= 1.20.0
- scipy >= 1.7.0
- pymatgen >= 2022.0.0
- spglib >= 1.16.0

## 注意事项

1. 确保已安装MOM2SSG包，因为POS2SSG依赖它
2. 输入POSCAR文件格式需要符合VASP标准
3. 磁原子类型和磁矩列表长度必须匹配
4. 建议使用`--verbose`选项查看详细处理过程

## 故障排除

如果遇到导入错误，请检查：
1. 是否正确安装了所有依赖项
2. MOM2SSG包是否可用
3. Python环境是否正确激活

如果遇到构建错误，请检查：
1. Python版本是否符合要求
2. 是否安装了构建工具（setuptools, wheel, build）
3. 是否有足够的磁盘空间
