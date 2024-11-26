<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Causal-DRF: Conditional Kernel Treatment Effect Estimation</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f9f9f9;
        }
        h1 {
            color: #333;
        }
        a {
            color: #1a73e8;
            text-decoration: none;
        }
        a:hover {
            text-decoration: underline;
        }
        code {
            background: #f4f4f4;
            padding: 2px 4px;
            border-radius: 4px;
            font-size: 0.95em;
        }
        pre {
            background: #f4f4f4;
            padding: 10px;
            border-radius: 4px;
            overflow-x: auto;
        }
    </style>
</head>
<body>
    <h1>Causal-DRF: Conditional Kernel Treatment Effect Estimation using Distributional Random Forest</h1>
    <p>
        This repository contains the code to replicate the results presented in the paper: 
        <strong>
            <a href="https://arxiv.org/abs/2411.08778" target="_blank">
                Causal-DRF: Conditional Kernel Treatment Effect Estimation using Distributional Random Forest
            </a>
        </strong>.
    </p>

    <h2>Requirements</h2>
    <p>
        To run this code, you will need to install a new version of the <code>drf</code> R package from the 
        <strong>causal-clean</strong> branch on GitHub.
    </p>

    <h2>Installation</h2>
    <p>Follow these steps to install the required version of the <code>drf</code> package:</p>
    <ol>
        <li>
            Clone the <code>drf</code> repository from the <strong>causal-clean</strong> branch:
            <pre><code>git clone --branch causal-clean https://github.com/herbps10/drf.git</code></pre>
        </li>
        <li>
            Navigate to the <code>r-package</code> directory:
            <pre><code>cd drf/r-package</code></pre>
        </li>
        <li>
            Build and install the R package:
            <pre><code>Rscript build_package.R</code></pre>
        </li>
    </ol>

    <h2>Usage</h2>
    <p>
        Once the <code>drf</code> package is installed, you can use the provided code to replicate the results 
        from the paper. Additional details about running the code and reproducing experiments will be provided 
        in the repository's scripts.
    </p>

    <h2>Citation</h2>
    <p>
        If you find this work useful, please consider citing the associated paper:
    </p>
    <blockquote>
        "Causal-DRF: Conditional Kernel Treatment Effect Estimation using Distributional Random Forest" <br>
        Available at: 
        <a href="https://arxiv.org/abs/2411.08778" target="_blank">https://arxiv.org/abs/2411.08778</a>
    </blockquote>

    <hr>
    <p>For further questions or issues, feel free to open an issue on this repository.</p>
</body>
</html>
![image](https://github.com/user-attachments/assets/f48495dc-b41c-46f5-8c9f-f8a502a184a9)
