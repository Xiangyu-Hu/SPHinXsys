# A5 shared Reload（不入库）

`L1|L2|L3/Reload.xml` 由 **relax** 生成，用于三端对齐初值；**不要提交 Git**。

生成并发布示例：

```bash
env RUN_SUFFIX=share01 GPU_PAPER_PUBLISH_SHARED_RELOAD=1 \
  bash docs/gpu_paper_discussion/scripts/run_a5.sh L2 gpu_fp32 relax
```

缺文件时，`run_a5.sh` 在 smoke/validate 前也会自动先 relax。
详见同级 `../NOTES.md`。
