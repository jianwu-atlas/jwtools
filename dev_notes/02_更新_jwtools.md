# jwtools 本地更新脚本 — `update_jwtools.sh`

> **用途**：一键重新构建 + 安装本地 jwtools 包，替代手动敲5行命令。
> **触发场景**：每次修改/新增 `jwtools/R/*.R` 中的函数后，必须运行此脚本，
> 才能让 `library(jwtools)` 加载到最新版本（详见下方"核心原理"）。

---

## 脚本位置

```
~/bin/update_jwtools.sh
```

## 脚本内容（最新可用版本）

```bash
#!/bin/bash
# ============================================================
# update_jwtools.sh — 一键重新构建+安装 jwtools
# 修复记录: upgrade参数改用FALSE
#          （新版remotes包要求布尔值，不再接受"never"字符串）
# ============================================================
JWTOOLS_PATH="/Users/Jack/Documents/GitHub/jwtools"

Rscript -e "
if ('jwtools' %in% loadedNamespaces()) unloadNamespace('jwtools')
devtools::document('${JWTOOLS_PATH}')
devtools::install('${JWTOOLS_PATH}', quick = FALSE, upgrade = FALSE)
cat('\n\u2705 jwtools \u5df2\u66f4\u65b0\u81f3\u7248\u672c:', as.character(packageVersion('jwtools')), '\n')
"
```

## 使用方法

```bash
chmod +x ~/bin/update_jwtools.sh   # 首次使用需授权（仅需一次）
~/bin/update_jwtools.sh             # 每次改完函数后运行这一句
```

---

## 曾踩过的坑（Troubleshooting Log）

### 坑1：`upgrade` 参数报错

**报错信息**：

```
! `upgrade` must be a single TRUE, FALSE, or NA
Backtrace:
 1. devtools::install(..., upgrade = "never")
 2.   cli::cli_abort("{.arg upgrade} must be a single TRUE, FALSE, or NA")
```

**原因**：`remotes`/`devtools` 包版本升级后，`upgrade` 参数不再接受字符串
（旧版本可用 `"never"`/`"always"`/`"ask"`），新版本只接受布尔值。

**修复**：

```r
# 旧写法（新版本报错）
devtools::install(path, upgrade = "never")

# 新写法（正确）
devtools::install(path, upgrade = FALSE)
```

**判断是否命中此坑**：如果 `devtools::document()` 那一步显示成功
（"Updating jwtools documentation"），但卡在 `install()` 报 `upgrade` 相关
错误，就是这个问题，直接改参数即可，与代码逻辑无关。

---

### 坑2：改了本地源码，`library(jwtools)` 加载的还是旧版本

**核心原理（务必记住）**：

```
本地源码 (jwtools/R/*.R)
     |
     |--- git push -----------> GitHub远程仓库
     |                          （只是代码备份，不影响本地R环境！）
     |
     +--- devtools::install() -> 本地R系统库 (.libPaths())
                                     |
                          library(jwtools) 加载的正是这里
```

**结论**：`git push` 到 GitHub 不会自动更新本地系统里已安装的R包。
必须手动跑 `update_jwtools.sh`（或 `devtools::install()`）才会生效。

**建议排查顺序**：

1. 改完函数 → 先跑 `update_jwtools.sh`
2. 重启R session（RStudio: `Cmd/Ctrl + Shift + F10`），避免旧版本缓存在内存里
3. `library(jwtools)` 后用 `packageVersion("jwtools")` 确认版本号
4. `exists("函数名")` 确认新函数确实存在

---

## 验证更新是否成功（标准检查清单）

```r
# 1. 重启R session后执行
library(jwtools)

# 2. 检查版本号（对应 DESCRIPTION 文件里的 Version 字段）
packageVersion("jwtools")

# 3. 检查新函数是否存在
exists("ct_proportion_analysis")   # 应返回 TRUE

# 4. 检查帮助文档是否同步更新
?ct_proportion_analysis
```

---

## 相关联的开发实践建议

每次改动后建议同步递增版本号：

```r
usethis::use_version("patch")   # 小改动: 0.1.0 -> 0.1.1
usethis::use_version("minor")   # 大功能: 0.1.0 -> 0.2.0
```

配合 git commit，保留清晰的变更记录：

```bash
cd /Users/Jack/Documents/GitHub/jwtools
git add R/ man/ NAMESPACE DESCRIPTION
git commit -m "feat: fix upgrade param in update_jwtools.sh (upgrade=FALSE)"
git push
```

---

## 相关文件

- 脚本本体: `~/bin/update_jwtools.sh`
- 包源码根目录: `/Users/Jack/Documents/GitHub/jwtools/`
- 索引文件: `jwtools/dev_notes/00_INDEX.md`（需手动补充本文件条目，
  或运行 `update_log_index.sh` 自动刷新）

---

_记录日期: 2025-11-20_
_记录人: Jian Wu_