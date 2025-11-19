<img src="https://github.com/shenyulu/easyclimate/blob/main/docs/source/_static/easyclimate-logo.svg?raw=true" alt="easyclimate">

<h2 align="center">一行のコードで気候データを分析する</h2>

<p align="center">
<a href="https://easyclimate.readthedocs.io/en/latest/"><strong>ドキュメント</strong> 「最新」</a> •
<a href="https://easyclimate.readthedocs.io/en/main/"><strong>ドキュメント</strong> 「メインブランチ」</a> •
<a href="https://shenyulu.github.io/easyclimate/"><strong>ドキュメント</strong> 「開発ブランチ」</a> •
<a href="https://github.com/shenyulu/easyclimate/blob/main/CONTRIBUTING.md"><strong>貢献</strong></a>
</p>

![PyPI - Version](https://img.shields.io/pypi/v/easyclimate)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/easyclimate)
![PyPI - Downloads](https://img.shields.io/pypi/dm/easyclimate)
[![codecov](https://codecov.io/gh/shenyulu/easyclimate/graph/badge.svg?token=CBG3IO5A5A)](https://codecov.io/gh/shenyulu/easyclimate)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/shenyulu/easyclimate/main.svg)](https://results.pre-commit.ci/latest/github/shenyulu/easyclimate/main)
[![Documentation Status](https://readthedocs.org/projects/easyclimate/badge/?version=latest)](https://easyclimate.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/465206111.svg)](https://zenodo.org/doi/10.5281/zenodo.10279567)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/shenyulu/easyclimate/main?labpath=docs%2Fexample)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:53294fc5285255eae60d8cf5d50ec86f28bb970a/)](https://archive.softwareheritage.org/swh:1:dir:53294fc5285255eae60d8cf5d50ec86f28bb970a;origin=https://doi.org/10.5281/zenodo.10279567;visit=swh:1:snp:85515406366d3bdbb39344f4fd07373a49cd1f42;anchor=swh:1:rel:12cab9c3f8cd0bd24229e4d9a5c778fe9b10df63;path=/shenyulu-easyclimate-b0b85e0/)

<div align="center">
<center><a href = "../README.md">English</a> / <a href = "README.zh_CN.md">简体中文</a> / 日本語</center>
</div>

## 👋 概要

**Easy Climate** は、気候学を一行のコードで分析することを目的とした Python パッケージです。

現在まで、Python >= 3.10 に対応しています。

### ✨ プロジェクトの目的

* 冗長なデータおよびグラフィック処理コードを排除する。
* 最適な計算速度を得るために、numpy と xarray の構文を使用する。
* Dask を使用してスケーラブルな並列処理をサポートし、TB または PB 単位の気候データを処理できる。
* オープンソースソフトウェアで、アプリケーションの柔軟性を向上させる。

### 🚀 インストール方法

`easyclimate` は [pip](https://pypi.org/project/pip/) を使ってインストールできます。

```bash
pip install easyclimate
```

インストールに関する詳細は、ドキュメントを見てください。

## 🤖 AI アシスタント（実験的）

[DeepWiki](https://deepwiki.com/shenyulu/easyclimate) は、Devin チームによって開発されたインテリジェントなAI 駆動型の GitHub リポジトリアシスタントです。ユーザーは自然言語処理を利用してコードリポジトリに直接質問することができ、詳細かつドキュメントレベルの回答を得ることができます。DeepWiki は複雑な問題の詳細な分析のためのディープリサーチモードをサポートしています。

このプロジェクトで DeepWiki を試してみる：

🔗 https://deepwiki.com/shenyulu/easyclimate

## 💫 参加方法

👩🏾‍💻 **プロジェクト開発への貢献：**
私たちの [貢献ガイド](https://github.com/shenyulu/easyclimate/blob/main/CONTRIBUTING.md) をお読みいただき、どのように手伝ったりフィードバックを提供したりできるかをご確認ください。

🧑🏾‍🤝‍🧑🏼 **行動規範：**
このプロジェクトには [行動規範](https://github.com/shenyulu/easyclimate/blob/main/CODE_OF_CONDUCT.md) が定められています。
このプロジェクトに参加することで、その条件に従うことに同意したことになります。

❤️ **プロジェクトのスポンサーになる**
貴重なフィードバックや不具合報告を通じてソフトウェアの改善に貢献してくださった皆様に、心より感謝申し上げます。このプロジェクトはコミュニティ主導で運営されており、スポンサーの皆様のご支援に深く感謝しております。スポンサー一覧は[こちら❤️](https://easyclimate.readthedocs.io/en/latest/sponsor.html)でご確認いただけます。

> **インポスター症候群に関する免責事項：**
> あなたの助けが必要です。**本当に。** 頭の中に「まだ準備ができていない」「貢献するにはスキルが足りない」と言う小さな声が聞こえるかもしれませんが、その声は間違っています。最も重要なのは、**コードを書くこと以外にも貢献する価値のある方法がたくさんある** ということです。

## 🤗 貢献者

たくさん貢献者に感謝します！

[![Contributors](https://contrib.rocks/image?repo=shenyulu/easyclimate)](https://github.com/shenyulu/easyclimate/graphs/contributors)

## 🪪 オープンソースライセンス

これは無料のソフトウェアです。**BSD 3 Clause License**の条件で再配布および変更が可能です。
このライセンスのコピーは [`LICENSE`](https://github.com/shenyulu/easyclimate/blob/main/LICENSE) に記載されています。

## 💎 スター履歴

[![Star History Chart](https://api.star-history.com/svg?repos=shenyulu/easyclimate&type=Date)](https://star-history.com/#shenyulu/easyclimate&Date)
