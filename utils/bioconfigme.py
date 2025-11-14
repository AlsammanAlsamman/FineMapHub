"""
bioconfigme - Configuration management for FineMapHub pipeline

Centralized access to analysis and software configuration.
All configuration access should go through this module.
"""

import yaml
import os
from pathlib import Path
from typing import Any, Dict, Optional, Sequence


class ConfigError(Exception):
    """Raised when configuration is invalid or missing."""
    pass


def _load_config_file(file_path: str) -> Dict[str, Any]:
    """Load a single YAML configuration file."""
    path = Path(file_path)
    if not path.exists():
        raise ConfigError(f"Configuration file not found: {file_path}")
    
    try:
        with open(path, 'r') as f:
            return yaml.safe_load(f) or {}
    except yaml.YAMLError as e:
        raise ConfigError(f"Failed to parse YAML file {file_path}: {e}")
    except Exception as e:
        raise ConfigError(f"Failed to read configuration file {file_path}: {e}")


def _get_config() -> Dict[str, Any]:
    """Load and merge all configuration files."""
    config_dir = Path("configs")
    
    # Load analysis and software configs
    analysis_config = _load_config_file(config_dir / "analysis.yml")
    software_config = _load_config_file(config_dir / "software.yml")
    
    return {
        "analysis": analysis_config,
        "software": software_config
    }


def get_results_dir() -> str:
    """Get the results directory path from analysis configuration."""
    config = _get_config()
    results_dir = config["analysis"].get("results_dir")
    
    if not results_dir:
        raise ConfigError("results_dir not found in analysis.yml")
    
    return str(results_dir)


def get_software_module(tool_name: str) -> str:
    """Get the module string for a software tool."""
    config = _get_config()
    software = config["software"]
    
    if tool_name not in software:
        raise ConfigError(f"Software tool '{tool_name}' not found in software.yml")
    
    tool_config = software[tool_name]
    module = tool_config.get("module")
    
    if not module:
        raise ConfigError(f"Module not defined for tool '{tool_name}' in software.yml")
    
    return str(module)


def get_analysis_value(path: Sequence[str], *, required: bool = True, default: Optional[Any] = None) -> Any:
    """Get a value from analysis configuration using dot notation path."""
    config = _get_config()
    
    # Navigate through the path
    value = config["analysis"]
    path_str = ".".join(path)
    
    try:
        for key in path:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                if required:
                    raise ConfigError(f"Configuration key '{path_str}' not found in analysis.yml")
                return default
    except (TypeError, KeyError):
        if required:
            raise ConfigError(f"Configuration key '{path_str}' not found in analysis.yml")
        return default
    
    return value


def get_software_params(tool_name: str) -> Dict[str, Any]:
    """Get parameters for a software tool."""
    config = _get_config()
    software = config["software"]
    
    if tool_name not in software:
        raise ConfigError(f"Software tool '{tool_name}' not found in software.yml")
    
    tool_config = software[tool_name]
    return tool_config.get("params", {})


def get_default_resources() -> Dict[str, Any]:
    """Get default resource requirements."""
    try:
        return get_analysis_value(["default_resources"], required=False, default={
            "mem_mb": 32000,
            "cores": 2,
            "time": "00:30:00"
        })
    except ConfigError:
        return {
            "mem_mb": 32000,
            "cores": 2,
            "time": "00:30:00"
        }


def get_gwas_table_name(target_analysis_name: str) -> str:
    """Get the GWAS table name for a target analysis."""
    return get_analysis_value(["target_analysis", target_analysis_name, "gwas_table"])


def get_gwas_table_info(gwas_table_name: str, info_key: str) -> Any:
    """Get information from a GWAS table by name and key (e.g., sample_size, case, control)."""
    gwas_tables = get_analysis_value(["gwastables"])
    
    for table in gwas_tables:
        if table.get("name") == gwas_table_name:
            if info_key in table:
                return table[info_key]
            else:
                raise ConfigError(f"Key '{info_key}' not found in GWAS table '{gwas_table_name}'")
    
    raise ConfigError(f"GWAS table '{gwas_table_name}' not found in configuration")


def get_gwas_sample_size(target_analysis_name: str) -> int:
    """Get the sample size for a target analysis by looking up its GWAS table."""
    gwas_table_name = get_gwas_table_name(target_analysis_name)
    return int(get_gwas_table_info(gwas_table_name, "sample_size"))