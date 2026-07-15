# jwtools

个人科研常用 R 小工具函数集合。目的很简单：以后不用每次开新脚本都从
以前的项目里翻出 `qs_save_workspace()` / `qs_load_workspace()` 之类的
函数复制粘贴 —— 装一次包，`library(jwtools)` 就能直接用，以后写了新的
小工具函数也往这里加。

## 现在包含什么

| 函数 | 作用 | 对应你原来的写法 |
|---|---|---|
| `qs_save_workspace()` | 把当前环境所有变量存成一个 `.qs2` 文件 | `save(list = ls(all=TRUE), file = "xxx.RData")` |
| `qs_load_workspace()` | 从 `.qs2` 文件恢复所有变量 | `load("xxx.RData")` |

比 `.RData` 更快、文件更小，尤其是环境里有 Seurat 对象这种大东西的时候。

## 安装

### 方式一：本地安装（现在就能用，不用等发布到 GitHub）

把这个 `jwtools/` 文件夹放到你服务器/电脑上任意位置，然后：

```r
install.packages("devtools")   # 如果还没装过
devtools::install_local("path/to/jwtools")
```

### 方式二：放到 GitHub 后远程安装（推荐，方便以后在服务器和 Mac 上同步更新）

1. 在 GitHub 建一个新仓库，比如 `jwtools`
2. 把这个文件夹的内容整个传上去（结构保持不变：`DESCRIPTION`、`NAMESPACE`、
   `R/`、`man/` 等都在仓库根目录）
3. 以后在任何一台机器上：

```r
install.packages("remotes")   # 如果还没装过
remotes::install_github("JackNg88/jwtools")
```

改了函数、推送到 GitHub 之后，其他机器上重新跑一次
`remotes::install_github("JackNg88/jwtools")` 就能更新到最新版。

## 使用

```r
library(jwtools)

# 保存整个工作环境（等价于原来的 save(list = ls(all=TRUE), file = ...)）
qs_save_workspace("core_WT_YMO_workspace.qs2", nthreads = 14)

# 如果某个大对象（比如 immune.combined）已经单独存过，可以排除掉避免重复占磁盘
qs_save_workspace(
  "core_WT_YMO_workspace.qs2",
  nthreads = 14,
  exclude = c("immune.combined")
)

# 或者按正则排除一批变量（比如所有 ggplot 对象 p_xxx、临时对象 tmp_xxx）
qs_save_workspace(
  "core_WT_YMO_workspace.qs2",
  nthreads = 14,
  exclude_pattern = "^p_|^tmp_|^fig[0-9]"
)

# 读回来
qs_load_workspace("core_WT_YMO_workspace.qs2", nthreads = 14)
```

函数文档：`?qs_save_workspace`、`?qs_load_workspace`

## 以后怎么加新函数（保持包结构一致）

1. 在 `R/` 目录下新建一个 `.R` 文件（比如 `R/slim_and_save.R`），把函数写进去，
   函数上方按同样格式加 roxygen2 注释（`#'` 开头的那些行，写清楚
   `@param`、`@return`、`@export`）。
2. 装好 devtools 后跑一次：

```r
devtools::document()   # 根据 roxygen 注释自动生成/更新 NAMESPACE 和 man/*.Rd
devtools::load_all()   # 本地重新加载测试
devtools::check()      # 可选，检查包有没有问题
```

3. 提交到 GitHub，其他机器上 `remotes::install_github()` 更新。

这样 `NAMESPACE` 和帮助文档就不用手写了，`devtools::document()` 全自动生成。

## 依赖

- [qs2](https://cran.r-project.org/package=qs2)（保存/读取用）
- 跑测试的话还需要 `testthat (>= 3.0.0)`（可选）

```r
install.packages("qs2")
```
